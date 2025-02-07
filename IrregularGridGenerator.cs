using System.Collections;
using System.Collections.Generic;
using UnityEngine;

#if UNITY_EDITOR
using UnityEditor;
#endif

/// <summary>
/// Attach to an empty GameObject. In Play Mode, click the Step 1..5 buttons (and Done at Step 6)
/// to generate:
///   1) Hex Points (using axial coords + ring groups)
///   2) Triangulate them ring-by-ring and sextant-by-sextant
///   3) Randomly merge triangles into quads
///   4) Subdivide each face
///   5) Relax the final mesh
/// 
/// This version ensures a proper hex shape by starting each ring at (q=-i, r=0) and walking
/// around 6 directions. We also avoid tuple deconstruction to fix older C# version issues.
/// </summary>
public class IrregularQuadGridGenerator : MonoBehaviour {
    public enum GenerationStep {
        None,
        HexPoints,       // Step 1
        Triangulate,     // Step 2
        RandomMerging,   // Step 3
        Subdivide,       // Step 4
        Relax,           // Step 5
        Completed
    }

    [Header("Hex Grid Settings")]
    [Tooltip("Number of rings around center. e.g. 5..61")]
    public int ringCount = 7;

    [Header("Relaxation Settings")]
    public int relaxationIterations = 100; // set to 100 iterations as desired

    [Header("Gizmos Display Settings")]
    public float gizmoPointRadius = 0.06f;

    public GenerationStep currentStep = GenerationStep.None;

    // ------------------ Internal Data ------------------
    private List<Vector3> allPoints = new List<Vector3>();

    // ringIndexLists[i] = list of point indices for ring i
    private List<List<int>> ringIndexLists = new List<List<int>>();

    // Our faces (each is a list of indices)
    private List<Face> faces = new List<Face>();

    // For tracking shared edges
    private Dictionary<Edge, List<int>> edgeToFaceMap = new Dictionary<Edge, List<int>>();

    // Points that are on the edge of the mesh
    private List<int> boundaryPoints = new List<int>();

    // Helper classes
    private class Face {
        public List<int> indices;
        public Face(List<int> idx) {
            indices = idx;
        }
    }

    private struct Edge {
        public int a;
        public int b;
        public Edge(int i1, int i2) {
            if (i1 < i2) {
                a = i1;
                b = i2;
            } else {
                a = i2;
                b = i1;
            }
        }
    }

    private void OnDrawGizmos() {
        if (!Application.isPlaying) return;

        // Draw points (red)
        Gizmos.color = Color.red;
        for (int i = 0; i < allPoints.Count; i++) {
            Vector3 p = allPoints[i];
            Gizmos.DrawSphere(p, gizmoPointRadius);
            /*
            #if UNITY_EDITOR
            // If you want labels above each sphere:
            Vector3 labelPos = p + Vector3.up * 0.2f;
            Handles.Label(labelPos, $"ID: {i}");
            #endif
            */
        }

        // Draw face edges (white)
        Gizmos.color = Color.white;
        foreach (var face in faces) {
            if (face == null) continue;
            int count = face.indices.Count;
            for (int i = 0; i < count; i++) {
                Vector3 vA = allPoints[face.indices[i]];
                Vector3 vB = allPoints[face.indices[(i + 1) % count]];
                Gizmos.DrawLine(vA, vB);
            }
        }
    }

    // -----------------------------------------------------
    // Inspector Buttons
    // -----------------------------------------------------
    [ContextMenu("Clear")]
    public void ClearAll() {
        allPoints.Clear();
        faces.Clear();
        edgeToFaceMap.Clear();
        ringIndexLists.Clear();
        boundaryPoints.Clear();
        currentStep = GenerationStep.None;
    }

    [ContextMenu("Step 1 - Generate Hex Points")]
    public void Step1GenerateHexPoints() {
        ClearAll();
        GenerateHexPoints(ringCount);
        currentStep = GenerationStep.HexPoints;
    }

    [ContextMenu("Step 2 - Triangulate (Sextant)")]
    public void Step2Triangulate() {
        if (currentStep != GenerationStep.HexPoints) return;
        TriangulateSextantWise();
        currentStep = GenerationStep.Triangulate;
    }

    [ContextMenu("Step 3 - Random Merging into Quads")]
    public void Step3RandomMerge() {
        if (currentStep != GenerationStep.Triangulate) return;
        MergeRandomTriangles();
        currentStep = GenerationStep.RandomMerging;
    }

    [ContextMenu("Step 4 - Subdivide Faces")]
    public void Step4Subdivide() {
        if (currentStep != GenerationStep.RandomMerging) return;
        SubdivideFaces();
        currentStep = GenerationStep.Subdivide;
    }

    [ContextMenu("Step 5 - Relax Mesh")]
    public void Step5Relax() {
        if (currentStep != GenerationStep.Subdivide) return;
        StartCoroutine(RelaxMesh(relaxationIterations));
        currentStep = GenerationStep.Relax;
    }


    [ContextMenu("Done (Step 6)")]
    public void Done() {
        if (currentStep != GenerationStep.Relax) return;
        currentStep = GenerationStep.Completed;
    }

    // -----------------------------------------------------
    // STEP 1: Generate a Big Hex of Points, ring by ring
    // -----------------------------------------------------
    private void GenerateHexPoints(int rings) {
        // Axial coords -> index
        Dictionary<(int,int), int> axialToIndex = new Dictionary<(int,int), int>();

        allPoints = new List<Vector3>();
        ringIndexLists = new List<List<int>>();
        boundaryPoints = new List<int>();


        // ring 0 = single center at (0,0)
        allPoints.Add(Vector3.zero);
        ringIndexLists.Add(new List<int> { 0 });
        axialToIndex[(0,0)] = 0;
        int currentIndex = 1;

        // directions around the ring (pointy-top hex)
        (int dq,int dr)[] directions = {
            (1, -1),
            (1, 0),
            (0, 1),
            (-1, 1),
            (-1, 0),
            (0, -1)
        };

        for (int i = 1; i <= rings; i++) {
            List<(Vector2 pos, float angle, (int,int) axial)> ringPoints = 
                new List<(Vector2, float, (int,int))>();

            // Start corner for ring i
            int qCurr = -i;
            int rCurr = 0;

            // First corner
            Vector2 startPos = AxialToWorld2D(qCurr, rCurr);
            float startAngle = Mathf.Atan2(startPos.y, startPos.x);
            ringPoints.Add((startPos, startAngle, (qCurr, rCurr)));

            // walk around the ring
            foreach (var d in directions) {
                for (int step = 0; step < i; step++) {
                    // skip last step of last direction (wrap overlap)
                    if (step == i - 1 && d == directions[5]) {
                        continue;
                    }
                    qCurr += d.dq;
                    rCurr += d.dr;

                    Vector2 pos = AxialToWorld2D(qCurr, rCurr);
                    float ang = Mathf.Atan2(pos.y, pos.x);
                    ringPoints.Add((pos, ang, (qCurr, rCurr)));
                }
            }

            // create ring i (store final indices)
            List<int> ringList = new List<int>();
            for (int z = 0; z < ringPoints.Count; z++) {
                var rp = ringPoints[z];
                Vector3 pos3D = new Vector3(rp.pos.x, 0f, -rp.pos.y);
                allPoints.Add(pos3D);
                axialToIndex[rp.axial] = currentIndex;
                ringList.Add(currentIndex);

                // Mark boundary if outermost ring
                if (i == rings) {
                    boundaryPoints.Add(currentIndex);
                }

                currentIndex++;
            }
            ringIndexLists.Add(ringList);
        }
    }

    private Vector2 AxialToWorld2D(int q, int r) {
        float size = 1f;
        float x = size * Mathf.Sqrt(3f) * (q + (r * 0.5f));
        float y = size * (1.5f * r);
        return new Vector2(x, y);
    }

    // -----------------------------------------------------
    // STEP 2: Triangulate (Sextant by Sextant)
    // -----------------------------------------------------
    private void TriangulateSextantWise() {
        faces.Clear();
        edgeToFaceMap.Clear();

        for (int i = 1; i < ringIndexLists.Count; i++) {
            var prevRing = ringIndexLists[i - 1];
            var currRing = ringIndexLists[i];

            if (i == 1) {
                // ring 1 => 6 points around center => 6 triangles
                int centerIndex = ringIndexLists[0][0];
                for (int sext = 0; sext < 6; sext++) {
                    int a = currRing[sext];
                    int b = currRing[(sext + 1) % currRing.Count];
                    AddTriangle(centerIndex, a, b);
                }
            } else {
                // i > 1
                for (int sext = 0; sext < 6; sext++) {
                    List<int> segPrev = new List<int>();
                    for (int k = 0; k < i; k++) {
                        segPrev.Add(prevRing[(sext * (i - 1) + k) % prevRing.Count]);
                    }

                    List<int> segCurr = new List<int>();
                    for (int k = 0; k < i + 1; k++) {
                        segCurr.Add(currRing[(sext * i + k) % currRing.Count]);
                    }

                    TriangulateRingSegment(segPrev, segCurr);
                }
            }
        }
    }

    private void TriangulateRingSegment(List<int> prevSeg, List<int> currSeg) {
        int pCount = prevSeg.Count;
        int cCount = currSeg.Count;

        int pIndex = 0;
        int cIndex = 0;

        while (pIndex < pCount && cIndex < cCount - 1) {
            // tri 1: (pIndex, cIndex, cIndex+1)
            AddTriangle(prevSeg[pIndex], currSeg[cIndex], currSeg[cIndex + 1]);

            // tri 2: (pIndex, cIndex+1, pIndex+1)
            if (pIndex + 1 < pCount) {
                AddTriangle(prevSeg[pIndex], currSeg[cIndex + 1], prevSeg[pIndex + 1]);
            }
            cIndex++;
            if (pIndex + 1 < pCount) {
                pIndex++;
            }
        }
    }

    private void AddTriangle(int i0, int i1, int i2) {
        if (i0 == i1 || i1 == i2 || i2 == i0) {
            return; // degenerate
        }
        Face f = new Face(new List<int> { i0, i1, i2 });
        faces.Add(f);
        int faceIdx = faces.Count - 1;
        RegisterEdge(i0, i1, faceIdx);
        RegisterEdge(i1, i2, faceIdx);
        RegisterEdge(i0, i2, faceIdx);
    }

    private void RegisterEdge(int a, int b, int faceIndex) {
        Edge e = new Edge(a, b);
        if (!edgeToFaceMap.ContainsKey(e)) {
            edgeToFaceMap[e] = new List<int>();
        }
        edgeToFaceMap[e].Add(faceIndex);
    }

    // -----------------------------------------------------
    // STEP 3: Randomly merge triangles into quads
    // -----------------------------------------------------
    private void MergeRandomTriangles() {
        List<Edge> candidateEdges = new List<Edge>(edgeToFaceMap.Keys);
        System.Random rng = new System.Random();
        // Shuffle edges
        candidateEdges.Sort((a, b) => rng.Next(2) == 0 ? -1 : 1);

        HashSet<int> dirtyFaces = new HashSet<int>();

        foreach (var e in candidateEdges) {
            var faceList = edgeToFaceMap[e];
            if (faceList.Count != 2) continue;

            int fA = faceList[0];
            int fB = faceList[1];

            if (dirtyFaces.Contains(fA) || dirtyFaces.Contains(fB)) continue;

            var faceA = faces[fA];
            var faceB = faces[fB];

            if (faceA.indices.Count != 3 || faceB.indices.Count != 3) continue;

            // Merge the two triangles
            List<int> newFace = new List<int>();
            int sharedA = -1;
            int sharedB = -1;

            foreach (int idx in faceA.indices) {
                if (faceB.indices.Contains(idx)) {
                    if (sharedA == -1) {
                        sharedA = idx;
                    } else if (sharedB == -1) {
                        sharedB = idx;
                    }
                } else {
                    newFace.Add(idx);
                }
            }
            newFace.Add(sharedA);
            foreach (int idx in faceB.indices) {
                if (!newFace.Contains(idx) && idx != sharedB) {
                    newFace.Add(idx);
                }
            }
            newFace.Add(sharedB);

            Face newFaceObj = new Face(newFace);
            faces.Add(newFaceObj);
            int newFaceIdx = faces.Count - 1;

            dirtyFaces.Add(fA);
            dirtyFaces.Add(fB);

            // Register edges for new face
            RegisterEdge(newFace[0], newFace[1], newFaceIdx);
            RegisterEdge(newFace[1], newFace[2], newFaceIdx);
            RegisterEdge(newFace[2], newFace[3], newFaceIdx);
            RegisterEdge(newFace[3], newFace[0], newFaceIdx);
        }

        // Remove merged faces
        faces.RemoveAll(f => dirtyFaces.Contains(faces.IndexOf(f)));
        RebuildEdgeMap();
    }

    private void RebuildEdgeMap() {
        edgeToFaceMap.Clear();
        for (int i = 0; i < faces.Count; i++) {
            Face f = faces[i];
            if (f == null) continue;

            int c = f.indices.Count;
            for (int k = 0; k < c; k++) {
                int a = f.indices[k];
                int b = f.indices[(k + 1) % c];
                RegisterEdge(a, b, i);
            }
        }
    }

    // -----------------------------------------------------
    // STEP 4: Subdivide (tri or quad) into smaller quads
    // -----------------------------------------------------
    private void SubdivideFaces() {
        List<Face> newFaces = new List<Face>();

        Dictionary<(float, float), int> coordToPointIndex = new Dictionary<(float, float), int>();




        for (int i = 0; i < faces.Count; i++) {
            Face f = faces[i];
            if (f == null) continue;

            int vCount = f.indices.Count;

            // Centroid
            Vector3 centroid = Vector3.zero;
            for (int k = 0; k < vCount; k++) {
                centroid += allPoints[f.indices[k]];
            }
            centroid /= vCount;
            int cIndex = allPoints.Count;
            allPoints.Add(centroid);

            // Mid-edge points
            List<int> midEdgeIndices = new List<int>();
            for (int e = 0; e < vCount; e++) {
                int iA = f.indices[e];
                int iB = f.indices[(e + 1) % vCount];
                Vector3 mid = 0.5f * (allPoints[iA] + allPoints[iB]);

                // check if the midpoint is already in the list
                if (coordToPointIndex.TryGetValue((mid.x, mid.z), out int  midIndex)) {
                    midEdgeIndices.Add(midIndex);
                } else {
                    midIndex = allPoints.Count;

                    allPoints.Add(mid);
                    midEdgeIndices.Add(midIndex);
                    coordToPointIndex[(mid.x, mid.z)] = midIndex;
                }


                // if both iA and iB are boundary, the midpoint is also boundary
                if (boundaryPoints.Contains(iA) && boundaryPoints.Contains(iB)) {
                    boundaryPoints.Add(midIndex);
                }
            }

            // Build new quads
            for (int e = 0; e < vCount; e++) {
                int cornerA = f.indices[e];
                int mA = midEdgeIndices[e];
                int mB = midEdgeIndices[(e - 1 + vCount) % vCount];

                Face newF = new Face(new List<int> {
                    cIndex,
                    mB,
                    cornerA,
                    mA
                });
                newFaces.Add(newF);
            }
        }
        faces = newFaces;
        RebuildEdgeMap();
    }

    // -----------------------------------------------------
    // STEP 5: Relax interior points (Square-Based Quad Relaxation)
    // -----------------------------------------------------
    private IEnumerator RelaxMesh(int iterations) {
        float relaxationSpeed = 0.07f;

        for (int iter = 0; iter < iterations; iter++) {
            Dictionary<int, Vector3> displacementAccum = new Dictionary<int, Vector3>();
            Dictionary<int, int> displacementCount = new Dictionary<int, int>();

            foreach (Face face in faces) {
                Vector3[] quadVerts = new Vector3[4];
                for (int i = 0; i < 4; i++) {
                    quadVerts[i] = allPoints[face.indices[i]];
                }

                Vector3 centroid = (quadVerts[0] + quadVerts[1] + quadVerts[2] + quadVerts[3]) / 4f;

                float bestDist = float.MaxValue;
                Vector3[] bestTarget = new Vector3[4];

                float area = CalculateQuadArea(quadVerts);
                float radius = Mathf.Sqrt(area)/Mathf.Sqrt(2);

                for (int i = 0; i < 4; i++) {
                    Vector3 direction = (quadVerts[i] - centroid).normalized * radius;
                    float dist = 0f;
                    for (int j = 0; j < 4; j++) {
                        int ind = (i + j) % 4;

                        // new target point
                        Vector3 newTarget = centroid + (Quaternion.AngleAxis(-90 * j, Vector3.up) * direction) * radius;
                        dist += (newTarget - quadVerts[ind]).magnitude;

                    }
                    if (dist < bestDist) {
                        bestDist = dist;

                        for (int j = 0; j < 4; j++) {
                            int ind = (i + j) % 4;
                            bestTarget[ind] = centroid + (Quaternion.AngleAxis(-90 * j, Vector3.up) * direction) * radius;
                        }
                    }
                }

                for (int i = 0; i < 4; i++) {
                    float dist = (bestTarget[i] - quadVerts[i]).magnitude;
                    Vector3 displacement = Vector3.Lerp(quadVerts[i], bestTarget[i], relaxationSpeed * Mathf.Abs((radius - dist))/radius) - quadVerts[i];
                    displacementAccum[face.indices[i]] = displacementAccum.GetValueOrDefault(face.indices[i], Vector3.zero) + displacement;
                    displacementCount[face.indices[i]] = displacementCount.GetValueOrDefault(face.indices[i], 0) + 1;
                }
            }


            int counter = 0;
            // Apply the displacement if not a boundary point
            foreach (var displacement in displacementAccum) {
                if (!boundaryPoints.Contains(displacement.Key)) {
                    int ind = displacement.Key;
                    allPoints[ind] += displacementAccum[ind] / displacementCount[ind];
                }
                counter++;
            }

            Debug.Log("Relaxed iteration " + (iter + 1) + " / " + iterations);
            // wait 0.05 seconds
            yield return new WaitForSeconds(0.05f);
        }
    }

    float CalculateQuadArea(Vector3[] quadVerts) {
        float area1 = Vector3.Cross(quadVerts[1] - quadVerts[0], quadVerts[2] - quadVerts[0]).magnitude;
        float area2 = Vector3.Cross(quadVerts[2] - quadVerts[0], quadVerts[3] - quadVerts[0]).magnitude;
        return (area1 + area2) * 0.5f;
    }
}
