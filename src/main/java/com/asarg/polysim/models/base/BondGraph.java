package com.asarg.polysim.models.base;

import javafx.util.Pair;
import java.util.*;
/**
 * Created by Angel A. Cantu in 2018
 */
public class BondGraph {
    private ArrayList<ArrayList<Vertex>> preFaceVertices = new ArrayList<ArrayList<Vertex>>();
    private ArrayList<faceVertex> face_vertices = new ArrayList<>();
    private int totalAmountOfEdges = 0;
    private List<Vertex> vertices = new ArrayList<>();


    public class Vertex {
        private Tile thisTile;
        private String vertexLabel;
        private Coordinate location;
        private ArrayList<Border> itsBoundary = new ArrayList<Border>();
        private List<Edge> neighbors = new ArrayList<>();
        private Edge northNeighbor;
        private Edge southNeighbor;
        private Edge westNeighbor;
        private Edge eastNeighbor;
        private int mark = 0;
        private int pathMark = 0;
        private int edgeAmount = 0;
        private int boundaryFlag = 0;
        public void addEdge(Tile from, Tile to, Vertex b, int strength, String where) {
            totalAmountOfEdges++;
            neighbors.add(new Edge(from, to, b, strength, where));
            Edge e = new Edge(from, to, b, strength, where);
            edgeAmount = edgeAmount + 1;
            if (where == "N") {
                northNeighbor = e;
            } else if (where == "S") {
                southNeighbor = e;
            } else if (where == "W") {
                westNeighbor = e;
            } else if (where == "E") {
                eastNeighbor = e;
            }
        }
        public Coordinate getLocation() {
            return location;
        }
        public Tile getTile() {
            return thisTile;
        }
    }
    public class Edge {
        private Tile tileFrom;
        private Tile tileTo;
        private Vertex toVertex;
        private String direction = "X";
        private int bondStrength = 0;
        private int inside = -1;
        private int taken = 0;
        private Edge(Tile from, Tile to, Vertex v, int strength, String dir) {
            tileFrom = from;
            tileTo = to;
            toVertex = v;
            bondStrength = strength;
            direction = dir;
        }
        public Tile getTileFrom() {
            return tileFrom;
        }

        public Tile getTileTo() {
            return tileTo;
        }
    }
    private class faceVertex{
        private int inside = -1;
        private String label;
        private ArrayList<Pair<String, String>> edges = new ArrayList<>();
        private ArrayList<faceEdge> neighbors = new ArrayList<>();
        private ArrayList<Integer> strengths = new ArrayList<>();
        private faceVertex(String name){
            label = name;
        }
    }
    private class faceEdge {
        private faceVertex toVertex;
        private String source;
        private String original_source;
        private String original_to;
        private int bondStrength = 0;
        private faceEdge(String x, String y, String from, faceVertex b, int strength) {
            source = from;
            original_source = x;
            original_to = y;
            toVertex = new faceVertex(b.label);
            bondStrength = strength;
        }
    }
    private class Border{
        private String ofPath = "X";
        private Pair<String, String> inAndOutBorder;
        private Border(String path, String in, String out){
            ofPath = path;
            inAndOutBorder = new Pair<>(in, out);
        }
        private boolean Compare(String x, String y){
            Pair<String, String> c = new Pair<>(x, y);
            if (inAndOutBorder.equals(c)){
                return true;
            }
            else {
                return false;
            }
        }
    }
    private class bellmanGraph{
        private int amountOfV;
        private int amountOfE;
        private bellmanEdge edge[];

        private class bellmanEdge {
            private int src;
            private int dest;
            private int weight;
            private int flag;
            private int org_src;
            private int org_dest;
            private int id;
            private int redundantFlag;
            private bellmanEdge() {
                src = 0;
                dest = 0;
                weight = 0;
                flag = 0;
                org_dest = 0;
                org_src = 0;
                id = 0;
                redundantFlag = 0;
            }
        }
        private bellmanGraph(int v, int e) {
            amountOfV = v;
            amountOfE = e;
            edge = new bellmanEdge[amountOfE];
            for (int i = 0; i < amountOfE; ++i)
                edge[i] = new bellmanEdge();
        }
    }
    public List<Vertex> getListOfVertices(){
        return vertices;
    }
    public void addBondVertex(Tile tile, Coordinate pt) {
        Vertex v = new Vertex();
        v.vertexLabel = Integer.toString(vertices.size() - 1);
        v.thisTile = tile;
        v.location = pt;
        vertices.add(v);
    }
    public List<Edge> getMinCut(int tao) {
        boolean validAssembly = handleLeafs();
        if (!validAssembly){
            return new ArrayList<Edge>();
        }
        startGetFaces();
        insideWho();
        makeFaceVertices();
        setUpEdgesDualGraph();
        connectFaceEdges(false);
        connectFaceEdges(true);
        removeRedundantEdges();
        List<Edge> cutMe = minCut(tao);

        return cutMe;
    }
    private boolean handleLeafs(){
        int amountNoFaces = 0;
        while (true) {
            int amount = markLeafs();
            if (amount >= 2) {
                break;
            }
            else if (amountNoFaces == vertices.size() || amount == -1){
                System.out.print("INVALID ASSEMBLY || NO FACES");
                return false;
            }
            amountNoFaces++;
        }
        return true;
    }
    private void removeRedundantEdges(){
        for (faceVertex fv : face_vertices){
            for (faceEdge e : fv.neighbors){
                for (faceEdge ee: fv.neighbors){
                    if (ee.original_to == e.original_source && ee.original_source == e.original_to){
                        ee.original_source = "deleteMe";
                    }
                }
            }
        }
        for (faceVertex fv : face_vertices){
            Iterator<faceEdge> ee = fv.neighbors.iterator();
            while (ee.hasNext()){
                faceEdge e = ee.next();
                if (e.original_source == "deleteMe"){
                    ee.remove();
                }
            }
        }
    }
    private List<Edge> minCut(int tao){
        List<Integer> s_Belonging = new ArrayList<>();
        List<Integer> t_Belonging = new ArrayList<>();
        List<Integer> path = new ArrayList<>();
        for (Vertex s: vertices){
            findFaceOf(s_Belonging, s);
            for (Vertex t : vertices){
                if (t.vertexLabel == s.vertexLabel){
                    t_Belonging.clear();
                }
                else{
                    List<Integer> s_disjunction = new ArrayList<>();
                    s_disjunction.addAll(s_Belonging);
                    findFaceOf(t_Belonging, t);
                    s_disjunction.removeAll(t_Belonging);

                    List<Integer> t_disjunction = new ArrayList<>();
                    t_disjunction.addAll(s_Belonging);
                    findFaceOf(t_Belonging, t);
                    t_disjunction.removeAll(t_Belonging);

                    for (Integer s_vertex : s_Belonging){
                        for (Integer f_vertex : t_Belonging) {
                            if (!(f_vertex.equals(s_vertex))) {
                                if (s_vertex == -1){
                                    s_vertex = face_vertices.size() - 1;
                                }
                                if (f_vertex == -1){
                                    f_vertex = face_vertices.size() - 1;
                                }
                                ArrayList<Pair<Integer, Integer>> cutPath = new ArrayList<>();
                                int hasCut = bellmanFord(cutPath, s_vertex, path, f_vertex, tao);
                                if (hasCut == 1){
                                    ArrayList<Edge> cutMe = new ArrayList<>();
                                    for (Pair i : cutPath){
                                        for (Vertex v : vertices){
                                            for (Edge e : v.neighbors){
                                                if (i.getValue().equals(Integer.parseInt(e.toVertex.vertexLabel))
                                                        && i.getKey().equals(Integer.parseInt(v.vertexLabel))){
                                                    cutMe.add(e);
                                                }
                                            }
                                        }
                                    }
                                    return cutMe;
                                }
                                else{
                                    cutPath.clear();
                                }
                            }
                        }
                    }
                }
                t_Belonging.clear();
            }
            s_Belonging.clear();
        }
        return new ArrayList();
    }
    private void findFaceOf(List<Integer> s_t_Belonging, Vertex s_t){
        for (faceVertex fc : face_vertices){
            for (faceEdge fe : fc.neighbors){
                if (s_t.vertexLabel.equals(fe.original_source) || s_t.vertexLabel.equals(fe.original_to)){
                    s_t_Belonging.add(Integer.parseInt(fc.label));
                }
            }
        }
        Set<Integer> removeDup = new HashSet<>();
        removeDup.addAll(s_t_Belonging);
        s_t_Belonging.clear();
        s_t_Belonging.addAll(removeDup);
    }
    private void makeBellmanGraph(bellmanGraph graphOfBellmanFord){
        int k = 0;
        for (faceVertex fc : face_vertices){
            for (faceEdge fe : fc.neighbors){
                int neg_source = Integer.parseInt(fe.source);
                int neg_to = Integer.parseInt(fe.toVertex.label);

                if (fe.source == "-1"){
                    neg_source = face_vertices.size() - 1;
                }
                if (fe.toVertex.label == "-1"){
                    neg_to = face_vertices.size() - 1;
                }
                graphOfBellmanFord.edge[k].src = neg_source;
                graphOfBellmanFord.edge[k].dest = neg_to;
                graphOfBellmanFord.edge[k].org_src = Integer.parseInt(fe.original_source);
                graphOfBellmanFord.edge[k].org_dest = Integer.parseInt(fe.original_to);
                graphOfBellmanFord.edge[k].weight = fe.bondStrength;
                graphOfBellmanFord.edge[k].id = k;
                k++;
            }
        }
        for (bellmanGraph.bellmanEdge e : graphOfBellmanFord.edge){
            for (bellmanGraph.bellmanEdge ee : graphOfBellmanFord.edge){
                if (ee.redundantFlag == 0 && e.org_src == ee.org_dest && e.org_dest == ee.org_src){
                    ee.id = e.id;
                    ee.redundantFlag = 1;
                }
            }
        }

    }
    private int bellmanFord(ArrayList cutPath, int source, List<Integer> path, int dest, int tao){
        ArrayList<ArrayList<Pair<Integer, Integer>>> bellman = new ArrayList<ArrayList<Pair<Integer,Integer>>>();
        for (faceVertex _ : face_vertices){
            bellman.add(new ArrayList<Pair<Integer, Integer>>());
        }
        bellmanGraph graphOfBellmanFord = new bellmanGraph(face_vertices.size(), totalAmountOfEdges);
        makeBellmanGraph(graphOfBellmanFord);
        int amountVertices = graphOfBellmanFord.amountOfV;
        int amountEdges = graphOfBellmanFord.amountOfE;
        int dist[] = new int[amountVertices];
        int pred[] = new int[amountVertices];
        ArrayList<Pair<Integer, Integer>> cameFrom = new ArrayList<Pair<Integer, Integer>>();
        for (int i=0; i < amountVertices; ++i){
            dist[i] = Integer.MAX_VALUE;
            pred[i] = -1;
            cameFrom.add(new Pair<>(0, 0));
        }
        dist[source] = 0;
        pred[source] = -2;

        int firstFlag;
        for (int i = 1; i < amountVertices; ++i){
            for (int j = 0; j < amountEdges; j++){
                int u = graphOfBellmanFord.edge[j].src;
                int v = graphOfBellmanFord.edge[j].dest;
                int weight = graphOfBellmanFord.edge[j].weight;
                int id = graphOfBellmanFord.edge[j].id;
                if (dist[u] != Integer.MAX_VALUE &&
                        dist[u] + weight < dist[v] && v != source) {
                    Pair<Integer, Integer> p = new Pair<>(id, u);
                    firstFlag = 0;
                    ArrayList<Integer> values = new ArrayList<>();
                    for (Pair p_ : bellman.get(v)){
                        if (Integer.valueOf((int)p_.getKey()) == id){
                            firstFlag = 1;
                        }
                        values.add(Integer.valueOf((int)p_.getValue()));
                    }

                    values.add(u);
                    for (Integer val : values) {
                        for (Pair p_ : bellman.get(val)) {
                            if (Integer.valueOf((int)p_.getKey()) == id){
                                firstFlag = 1;
                            }
                        }
                    }
                    ArrayList<Pair<Integer, Integer>> holder = new ArrayList<Pair<Integer, Integer>>();
                    if (firstFlag == 0){
                        for (Pair p_ : bellman.get(v)){
                            if (!(Integer.valueOf((int)p_.getValue()) == u)){
                                holder.add(p_);
                            }
                        }
                        bellman.set(v, holder);
                        bellman.get(v).add(p);
                    }
                    if (firstFlag == 0){
                        dist[v] = dist[u] + weight;
                        pred[v] = u;
                    }

                    cameFrom.set(v, new Pair<>(graphOfBellmanFord.edge[j].org_dest, graphOfBellmanFord.edge[j].org_src));
                    graphOfBellmanFord.edge[j].flag = 2;
                    for (bellmanGraph.bellmanEdge e : graphOfBellmanFord.edge){
                        if (e.org_src == graphOfBellmanFord.edge[j].org_dest && e.org_dest == graphOfBellmanFord.edge[j].org_src){
                            e.flag = 2;
                        }
                    }
                    for (bellmanGraph.bellmanEdge e : graphOfBellmanFord.edge){
                        int destroyFlag = 0;
                        Pair<Integer, Integer> destroyer = new Pair<>(e.org_dest, e.org_src);
                        Pair<Integer, Integer> reyortsed = new Pair<>(e.org_src, e.org_dest);
                        for (Pair pp : cameFrom){
                            if (pp.getKey() == destroyer.getKey() && pp.getValue() == destroyer.getValue()
                                    || pp.getKey() == reyortsed.getKey() && pp.getValue() == reyortsed.getValue()){
                                destroyFlag = 1;
                            }
                        }
                        if (destroyFlag == 1){
                            e.flag = 2;
                        }
                        else{
                            e.flag = 0;
                        }
                    }
                }
            }
        }
        int iterator = 0;
        Iterator<ArrayList<Pair<Integer, Integer>>> bellIterator = bellman.iterator();
        ArrayList<Integer> bellBuffer = new ArrayList<>();
        while (bellIterator.hasNext()){
            ArrayList<Pair<Integer, Integer>> nextBellIt = bellIterator.next();
            Iterator<Pair<Integer, Integer>> bellItIterator = nextBellIt.iterator();
            while(bellItIterator.hasNext()){
                Pair<Integer, Integer> po = bellItIterator.next();
                if (po.getValue() == pred[iterator]){
                    bellBuffer.add(po.getKey());
                }
            }
            iterator++;
        }
        ArrayList<Integer> realPath = new ArrayList<>();
        for (bellmanGraph.bellmanEdge e : graphOfBellmanFord.edge){
            for (Integer i : bellBuffer){
                if (e.id == i){
                    e.flag = 1;
                    realPath.add(e.id);
                }
            }
        }
        dist[source] = Integer.MAX_VALUE;
        int cycleFlag = 0;
        for (int j = 0; j < amountEdges; ++j)
        {
            int u = graphOfBellmanFord.edge[j].src;
            int v = graphOfBellmanFord.edge[j].dest;
            int weight = graphOfBellmanFord.edge[j].weight;
            int f = graphOfBellmanFord.edge[j].flag;
            int id = graphOfBellmanFord.edge[j].id;
            if (u != source && dist[u] != Integer.MAX_VALUE &&
                    dist[u] + weight < tao  && f != 1 && v == source) {
                dist[v] = dist[u] + weight;
                pred[v] = u;
                Pair<Integer, Integer> pp = new Pair<>(id, u);
                bellman.get(v).add(pp);
                realPath.add(id);
                cycleFlag = 1;
                break;
            }
        }

        path.clear();
        int secondFlag = 0;
        for (Integer _ : pred){
            int previous = pred[dest];
            if (previous == Integer.MAX_VALUE){
                path.add(dest);
                if (dest == source){
                    secondFlag = 1;
                }
                break;
            }
            else{
                path.add(previous);
                dest = previous;
                if (dest == source){
                    secondFlag = 1;
                }
            }
            dest = previous;
            if (dest == -2){
                return 0;
            }
        }
        ArrayList<Integer> buffer = new ArrayList<>();
        ArrayList<Integer> idPath = new ArrayList<>();
        for (int j = 1; j < path.size(); j++){
            buffer.add(path.get(j));
        }
        for (int s = 0; s < path.size(); s++){
            if (s < buffer.size()){
                Iterator<Pair<Integer, Integer>> com = bellman.get(path.get(s)).iterator();
                while(com.hasNext()){
                    Pair<Integer, Integer> po = com.next();
                    if (po.getValue() == buffer.get(s)){
                        idPath.add(po.getKey());
                    }
                }
            }
            else{
                Iterator<Pair<Integer, Integer>> com = bellman.get(path.get(s)).iterator();
                while(com.hasNext()){
                    Pair<Integer, Integer> po = com.next();
                    if (po.getValue() == pred[path.get(s)]){
                        idPath.add(po.getKey());
                    }
                }
            }
        }
        Set<Integer> removeDup = new HashSet<>();
        removeDup.addAll(idPath);
        idPath.clear();
        idPath.addAll(removeDup);
        ArrayList<Pair<Integer, Integer>> finalPath = new ArrayList<>();
        for (Integer i : idPath){
            for (bellmanGraph.bellmanEdge e : graphOfBellmanFord.edge){
                if (i == e.id){
                    finalPath.add(new Pair<>(e.org_dest, e.org_src));
                }
            }
        }
        ArrayList<Integer> indices = new ArrayList<>();
        int firstCounter = 0;
        int flag = 0;
        for (Pair p : finalPath){
            int secondCounter = 0;
            for (Pair pp : finalPath){
                if (firstCounter == secondCounter){

                }
                else if (p.getValue() == pp.getKey() && p.getKey()== pp.getValue()){
                    for (Integer i : indices){
                        if (i == secondCounter){
                            flag = 1;
                        }
                    }
                    if (flag == 0){
                        indices.add(firstCounter);
                    }
                    flag = 0;
                }
                secondCounter++;
            }
            firstCounter++;
        }

        if (cycleFlag == 1 && dist[source] < tao && secondFlag != 0){
            for (Integer i : indices){
                cutPath.add(finalPath.get(i));
            }
            return 1;
        }
        else{
            return 0;
        }
    }
    private void connectFaceEdges(boolean x){
        if (!x){
            faceVertex outsideFace = new faceVertex("-1");
            face_vertices.add(outsideFace);
        }
        for (Vertex v : vertices){
            int flag = 0;
            if (v.northNeighbor != null && v.northNeighbor.direction != null
                    && !v.northNeighbor.direction.isEmpty()
                    && v.northNeighbor.taken == 0){
                Pair<String, String> check_pair = new Pair<>(v.vertexLabel, v.northNeighbor.toVertex.vertexLabel);
                Pair<String, String> r_check_pair = new Pair<>(v.northNeighbor.toVertex.vertexLabel, v.vertexLabel);
                for (faceVertex fc : face_vertices){
                    for (Pair fe : fc.edges){
                        if (check_pair.equals(fe)){
                            if (v.northNeighbor.inside == -1){
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.northNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(face_vertices.size() - 1), v.northNeighbor.bondStrength));
                                face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.northNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, fc, v.northNeighbor.bondStrength));
                            }
                            else{
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.northNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(v.northNeighbor.inside), v.northNeighbor.bondStrength));
                                face_vertices.get(v.northNeighbor.inside).neighbors.add(new faceEdge(v.northNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.northNeighbor.inside).label, fc, v.northNeighbor.bondStrength));
                            }
                        }
                    }
                }
                if (flag == 0){
                    for (faceVertex fc : face_vertices){
                        for (Pair fe : fc.edges){
                            if (r_check_pair.equals(fe)){
                                flag = 1;
                            }
                        }
                    }
                    if (flag == 0){
                        if (v.northNeighbor.inside == -1){
                            face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.northNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, face_vertices.get(face_vertices.size() - 1), v.northNeighbor.bondStrength));
                        }
                        else{
                            face_vertices.get(v.northNeighbor.inside).neighbors.add(new faceEdge(v.northNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.northNeighbor.inside).label, face_vertices.get(v.northNeighbor.inside), v.northNeighbor.bondStrength));
                        }

                    }
                }
                v.northNeighbor.taken = 1;
            }
            if (v.eastNeighbor != null && v.eastNeighbor.direction != null
                    && !v.eastNeighbor.direction.isEmpty()
                    && v.eastNeighbor.taken == 0){
                Pair<String, String> check_pair = new Pair<>(v.vertexLabel, v.eastNeighbor.toVertex.vertexLabel);
                Pair<String, String> r_check_pair = new Pair<>(v.eastNeighbor.toVertex.vertexLabel, v.vertexLabel);
                for (faceVertex fc : face_vertices){
                    for (Pair fe : fc.edges){
                        if (check_pair.equals(fe)){
                            if (v.eastNeighbor.inside == -1){
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.eastNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(face_vertices.size() - 1), v.eastNeighbor.bondStrength));
                                face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.eastNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, fc, v.eastNeighbor.bondStrength));
                            }
                            else{
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.eastNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(v.eastNeighbor.inside), v.eastNeighbor.bondStrength));
                                face_vertices.get(v.eastNeighbor.inside).neighbors.add(new faceEdge(v.eastNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.eastNeighbor.inside).label, fc, v.eastNeighbor.bondStrength));
                            }
                        }
                    }
                }
                if (flag == 0){
                    for (faceVertex fc : face_vertices){
                        for (Pair fe : fc.edges){
                            if (r_check_pair.equals(fe)){
                                flag = 1;
                            }
                        }
                    }
                    if (flag == 0){
                        if (v.eastNeighbor.inside == -1){
                            face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.eastNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, face_vertices.get(face_vertices.size() - 1), v.eastNeighbor.bondStrength));
                        }
                        else{
                            face_vertices.get(v.eastNeighbor.inside).neighbors.add(new faceEdge(v.eastNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.eastNeighbor.inside).label, face_vertices.get(v.eastNeighbor.inside), v.eastNeighbor.bondStrength));
                        }

                    }
                }
                v.eastNeighbor.taken = 1;
            }
            if (v.westNeighbor != null && v.westNeighbor.direction != null
                    && !v.westNeighbor.direction.isEmpty()
                    && v.westNeighbor.taken == 0){
                Pair<String, String> check_pair = new Pair<>(v.vertexLabel, v.westNeighbor.toVertex.vertexLabel);
                Pair<String, String> r_check_pair = new Pair<>(v.westNeighbor.toVertex.vertexLabel, v.vertexLabel);
                for (faceVertex fc : face_vertices){
                    for (Pair fe : fc.edges){
                        if (check_pair.equals(fe)){
                            if (v.westNeighbor.inside == -1){
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.westNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(face_vertices.size() - 1), v.westNeighbor.bondStrength));
                                face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.westNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, fc, v.westNeighbor.bondStrength));
                            }
                            else{
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.westNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(v.westNeighbor.inside), v.westNeighbor.bondStrength));
                                face_vertices.get(v.westNeighbor.inside).neighbors.add(new faceEdge(v.westNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.westNeighbor.inside).label, fc, v.westNeighbor.bondStrength));
                            }
                        }
                    }
                }
                if (flag == 0){
                    for (faceVertex fc : face_vertices){
                        for (Pair fe : fc.edges){
                            if (r_check_pair.equals(fe)){
                                flag = 1;
                            }
                        }
                    }
                    if (flag == 0){
                        if (v.westNeighbor.inside == -1){
                            face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.westNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, face_vertices.get(face_vertices.size() - 1), v.westNeighbor.bondStrength));
                        }
                        else{
                            face_vertices.get(v.westNeighbor.inside).neighbors.add(new faceEdge(v.westNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.westNeighbor.inside).label, face_vertices.get(v.westNeighbor.inside), v.westNeighbor.bondStrength));
                        }

                    }
                }
                v.westNeighbor.taken = 1;
            }
            if (v.southNeighbor != null && v.southNeighbor.direction != null
                    && !v.southNeighbor.direction.isEmpty()
                    && v.southNeighbor.taken == 0){
                Pair<String, String> check_pair = new Pair<>(v.vertexLabel, v.southNeighbor.toVertex.vertexLabel);
                Pair<String, String> r_check_pair = new Pair<>(v.southNeighbor.toVertex.vertexLabel, v.vertexLabel);
                for (faceVertex fc : face_vertices){
                    for (Pair fe : fc.edges){
                        if (check_pair.equals(fe)){
                            if (v.southNeighbor.inside == -1){
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.southNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(face_vertices.size() - 1), v.southNeighbor.bondStrength));
                                face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.southNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, fc, v.southNeighbor.bondStrength));
                            }
                            else{
                                flag = 1;
                                fc.neighbors.add(new faceEdge(v.vertexLabel, v.southNeighbor.toVertex.vertexLabel, fc.label, face_vertices.get(v.southNeighbor.inside), v.southNeighbor.bondStrength));
                                face_vertices.get(v.southNeighbor.inside).neighbors.add(new faceEdge(v.southNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.southNeighbor.inside).label, fc, v.southNeighbor.bondStrength));
                            }
                        }
                    }
                }
                if (flag == 0){
                    for (faceVertex fc : face_vertices){
                        for (Pair fe : fc.edges){
                            if (r_check_pair.equals(fe)){
                                flag = 1;
                            }
                        }
                    }
                    if (flag == 0){
                        if (v.southNeighbor.inside == -1){
                            face_vertices.get(face_vertices.size() - 1).neighbors.add(new faceEdge(v.southNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(face_vertices.size() - 1).label, face_vertices.get(face_vertices.size() - 1), v.southNeighbor.bondStrength));
                        }
                        else{
                            face_vertices.get(v.southNeighbor.inside).neighbors.add(new faceEdge(v.southNeighbor.toVertex.vertexLabel, v.vertexLabel, face_vertices.get(v.southNeighbor.inside).label, face_vertices.get(v.southNeighbor.inside), v.southNeighbor.bondStrength));
                        }

                    }
                }
                v.southNeighbor.taken = 1;
            }
        }
    }
    private void makeFaceVertices() {
        int i = 0;
        while(i < preFaceVertices.size()){
            faceVertex face_vertex = new faceVertex(Integer.toString(i));
            face_vertices.add(face_vertex);
            makeFaceEdges(i);
            i++;
        }
    }
    private void makeFaceEdges(int i){
        for (int j = 0; j < preFaceVertices.get(i).size() - 1; j++){

            if (preFaceVertices.get(i).get(j).northNeighbor != null && preFaceVertices.get(i).get(j).northNeighbor.direction != null
                    && !preFaceVertices.get(i).get(j).northNeighbor.direction.isEmpty()
                    && preFaceVertices.get(i).get(j).northNeighbor.inside != -1){
                if (face_vertices.get(i).inside < preFaceVertices.get(i).get(j).northNeighbor.inside){
                    face_vertices.get(i).inside = preFaceVertices.get(i).get(j).northNeighbor.inside;
                }
                face_vertices.get(i).strengths.add(preFaceVertices.get(i).get(j).northNeighbor.bondStrength);
            }
            if (preFaceVertices.get(i).get(j).eastNeighbor != null && preFaceVertices.get(i).get(j).eastNeighbor.direction != null
                    && !preFaceVertices.get(i).get(j).eastNeighbor.direction.isEmpty()
                    && preFaceVertices.get(i).get(j).eastNeighbor.inside != -1){
                if (face_vertices.get(i).inside < preFaceVertices.get(i).get(j).eastNeighbor.inside){
                    face_vertices.get(i).inside = preFaceVertices.get(i).get(j).eastNeighbor.inside;
                }
                face_vertices.get(i).strengths.add(preFaceVertices.get(i).get(j).eastNeighbor.bondStrength);
            }
            if (preFaceVertices.get(i).get(j).westNeighbor != null && preFaceVertices.get(i).get(j).westNeighbor.direction != null
                    && !preFaceVertices.get(i).get(j).westNeighbor.direction.isEmpty()
                    && preFaceVertices.get(i).get(j).westNeighbor.inside != -1){
                if (face_vertices.get(i).inside < preFaceVertices.get(i).get(j).westNeighbor.inside){
                    face_vertices.get(i).inside = preFaceVertices.get(i).get(j).westNeighbor.inside;
                }
                face_vertices.get(i).strengths.add(preFaceVertices.get(i).get(j).westNeighbor.bondStrength);
            }
            if (preFaceVertices.get(i).get(j).southNeighbor != null && preFaceVertices.get(i).get(j).southNeighbor.direction != null
                    && !preFaceVertices.get(i).get(j).southNeighbor.direction.isEmpty()
                    && preFaceVertices.get(i).get(j).southNeighbor.inside != -1){
                if (face_vertices.get(i).inside < preFaceVertices.get(i).get(j).southNeighbor.inside){
                    face_vertices.get(i).inside = preFaceVertices.get(i).get(j).southNeighbor.inside;
                }
                face_vertices.get(i).strengths.add(preFaceVertices.get(i).get(j).southNeighbor.bondStrength);
            }
            Pair<String, String> pair = new Pair(preFaceVertices.get(i).get(j).vertexLabel, preFaceVertices.get(i).get(j + 1).vertexLabel);
            face_vertices.get(i).edges.add(pair);
        }
    }
    private void setUpEdgesDualGraph(){
        int iterator = 0;
        for (faceVertex fv : face_vertices){
            for (Pair e : fv.edges){
                preConnectEdges(e, iterator);
            }
            iterator++;
        }
    }
    private void preConnectEdges(Pair e,  int which) {
        int iterator = 0;
        Pair<String, String> e_reversed = new Pair(e.getValue(), e.getKey());
        for (faceVertex fc : face_vertices) {
            if (!(iterator == which)) {
                for (Pair vv : fc.edges) {
                    if (e_reversed.equals(vv) || e.equals(vv)){
                        preCreateEdgeBetweenFaces(vv.getKey().toString(), vv.getValue().toString(), which, iterator);
                    }
                }
            }
            iterator++;
        }
    }
    private void preCreateEdgeBetweenFaces(String x, String y, int which, int count){
        int str = 0;
        String firstName = "";
        String secondName = "";
        for (Vertex v : vertices){
            if (v.vertexLabel == x) {
                if (v.northNeighbor != null && v.northNeighbor.direction != null
                        && !v.northNeighbor.direction.isEmpty()
                        && v.northNeighbor.toVertex.vertexLabel.equals(y)) {
                    str = v.northNeighbor.bondStrength;
                    v.northNeighbor.taken = 1;
                    firstName = v.vertexLabel;
                    secondName = v.northNeighbor.toVertex.vertexLabel;
                }
                if (v.eastNeighbor != null && v.eastNeighbor.direction != null
                        && !v.eastNeighbor.direction.isEmpty()
                        && v.eastNeighbor.toVertex.vertexLabel.equals(y)) {
                    str = v.eastNeighbor.bondStrength;
                    v.eastNeighbor.taken = 1;
                    firstName = v.vertexLabel;
                    secondName = v.eastNeighbor.toVertex.vertexLabel;
                }
                if (v.westNeighbor != null && v.westNeighbor.direction != null
                        && !v.westNeighbor.direction.isEmpty()
                        && v.westNeighbor.toVertex.vertexLabel.equals(y)) {
                    str = v.westNeighbor.bondStrength;
                    v.westNeighbor.taken = 1;
                    firstName = v.vertexLabel;
                    secondName = v.westNeighbor.toVertex.vertexLabel;
                }
                if (v.southNeighbor != null && v.southNeighbor.direction != null
                        && !v.southNeighbor.direction.isEmpty()
                        && v.southNeighbor.toVertex.vertexLabel.equals(y)) {
                    str = v.southNeighbor.bondStrength;
                    v.southNeighbor.taken = 1;
                    firstName = v.vertexLabel;
                    secondName = v.southNeighbor.toVertex.vertexLabel;
                }
                break;
            }
        }
        faceEdge e = new faceEdge(firstName, secondName, Integer.toString(which), face_vertices.get(count), str);
        face_vertices.get(which).neighbors.add(e);
    }
    private void markBoundary(List<Vertex> path, boolean which){
        if (which){
            for (Vertex v : path){
                v.boundaryFlag = 1;
            }
        }
        else{
            for (Vertex v : path){
                v.boundaryFlag = 0;
            }
        }
    }
    private int getIndexOfAttribute(Vertex v, int boundaryNum){
        if (!v.itsBoundary.isEmpty()){
            for (int i = 0; i < v.itsBoundary.size(); i++) {
                if (v.itsBoundary.get(i).ofPath.equals(Integer.toString(boundaryNum))) {
                    return i;
                }
            }
            return -1;
        }
        else{
            return 0;
        }
    }
    private void insideWho(){
        int boundaryNum = 0;
        int flag;
        int index;
        for (ArrayList<Vertex> v : preFaceVertices){
            markBoundary(v, true);
            for(Vertex vv : v){
                flag = 0;
                if (flag == 0 && vv.eastNeighbor != null && vv.eastNeighbor.direction != null
                        && !vv.eastNeighbor.direction.isEmpty() && vv.eastNeighbor.toVertex.boundaryFlag != 1){
                    index = getIndexOfAttribute(vv, boundaryNum);
                    if (index != -1){
                        if (vv.itsBoundary.get(index).Compare("E", "W")){
                            vv.eastNeighbor.inside = boundaryNum;
                            DFSmarking(vv.eastNeighbor.toVertex, boundaryNum);
                            flag = 1;
                        }
                    }
                }
                if (flag == 0 && vv.westNeighbor != null && vv.westNeighbor.direction != null
                        && !vv.westNeighbor.direction.isEmpty() && vv.westNeighbor.toVertex.boundaryFlag != 1){
                    index = getIndexOfAttribute(vv, boundaryNum);
                    if (index != -1){
                        if (vv.itsBoundary.get(index).Compare("W", "E")){
                            vv.westNeighbor.inside = boundaryNum;
                            DFSmarking(vv.westNeighbor.toVertex, boundaryNum);
                            flag = 1;
                        }
                    }
                }

                if (flag == 0 && vv.northNeighbor != null && vv.northNeighbor.direction != null
                        && !vv.northNeighbor.direction.isEmpty() && vv.northNeighbor.toVertex.boundaryFlag != 1){
                    index = getIndexOfAttribute(vv, boundaryNum);
                    if (index != -1){
                        if (vv.itsBoundary.get(index).Compare("N", "S")){
                            vv.northNeighbor.inside = boundaryNum;
                            DFSmarking(vv.northNeighbor.toVertex, boundaryNum);
                            flag = 1;
                        }
                    }
                }
                if (flag == 0 && vv.southNeighbor != null && vv.southNeighbor.direction != null
                        && !vv.southNeighbor.direction.isEmpty() && vv.southNeighbor.toVertex.boundaryFlag != 1){
                    index = getIndexOfAttribute(vv, boundaryNum);
                    if (index != -1){
                        if (vv.itsBoundary.get(index).Compare("S", "N")){
                            vv.southNeighbor.inside = boundaryNum;
                            DFSmarking(vv.southNeighbor.toVertex, boundaryNum);
                        }
                    }
                }
            }
            markBoundary(v, false);
            boundaryNum++;
        }
    }
    private void DFSmarking(Vertex v, int inside){
        ArrayList<Vertex> insiders = new ArrayList<>();
        insiders.add(v);
        do {
            for (Vertex vv : insiders) {
                if (vv.southNeighbor != null && vv.southNeighbor.direction != null
                        && !vv.southNeighbor.direction.isEmpty() && vv.southNeighbor.inside != inside) {
                    vv.southNeighbor.inside = inside;
                    if (vv.southNeighbor.toVertex.boundaryFlag != 1) {
                        insiders.add(vv.southNeighbor.toVertex);
                    }
                }
                if (vv.eastNeighbor != null && vv.eastNeighbor.direction != null
                        && !vv.eastNeighbor.direction.isEmpty() && vv.eastNeighbor.inside != inside) {
                    vv.eastNeighbor.inside = inside;
                    if (vv.eastNeighbor.toVertex.boundaryFlag != 1) {
                        insiders.add(vv.eastNeighbor.toVertex);
                    }
                }
                if (vv.westNeighbor != null && vv.westNeighbor.direction != null
                        && !vv.westNeighbor.direction.isEmpty() && vv.westNeighbor.inside != inside) {
                    vv.westNeighbor.inside = inside;
                    if (vv.westNeighbor.toVertex.boundaryFlag != 1) {
                        insiders.add(vv.westNeighbor.toVertex);
                    }
                }
                if (vv.northNeighbor != null && vv.northNeighbor.direction != null
                        && !vv.northNeighbor.direction.isEmpty() && vv.northNeighbor.inside != inside) {
                    vv.northNeighbor.inside = inside;
                    if (vv.northNeighbor.toVertex.boundaryFlag != 1) {
                        insiders.add(vv.northNeighbor.toVertex);
                    }
                }
                insiders.remove(insiders.indexOf(vv));
                break;
            }
        }while(!insiders.isEmpty());
    }
    private void startGetFaces() {
        ArrayList<Vertex> faces = new ArrayList<Vertex>();
        for (Vertex v : vertices) {
            if (!(v.neighbors.size() == 1)) {
                for (Edge e : v.neighbors) {
                    if (e.direction.equals("N")) {
                        faces = getFaces(v, e);
                        if (!faces.isEmpty()){
                            preFaceVertices.add(faces);
                        }
                    }
                }
            }
        }
        if (!preFaceVertices.isEmpty()){
            int flag = 0;
            ArrayList<ArrayList<Vertex>> path = new ArrayList<ArrayList<Vertex>>();
            path.add(preFaceVertices.get(0));
            for (ArrayList<Vertex> firstPath : preFaceVertices){
                for (ArrayList<Vertex> secondPath : path){
                    if (listEqualsIgnoreOrder(secondPath, firstPath)){
                        flag = 1;
                    }
                }
                if (flag == 0){
                    path.add(firstPath);
                }
                else{
                    flag = 0;
                }
            }
            preFaceVertices = path;
        }
        Collections.sort(preFaceVertices, new Comparator<ArrayList>(){
            public int compare(ArrayList a1, ArrayList a2) {
                return a2.size() - a1.size();
            }
        });
        buildBoundaries();
    }
    private void buildBoundaries(){
        int i = 0;
        int j = 0;
        for (ArrayList<Vertex> path : preFaceVertices){
            for (Vertex v : path){
                if (!(i >= path.size() - 1)){
                    if (v.northNeighbor != null && v.northNeighbor.direction != null
                            && !v.northNeighbor.direction.isEmpty() && v.northNeighbor.toVertex == path.get(i + 1)){
                        v.itsBoundary.add(new Border(Integer.toString(j), "W", "E"));
                    }
                    else if (v.westNeighbor != null && v.westNeighbor.direction != null
                            && !v.westNeighbor.direction.isEmpty() && v.westNeighbor.toVertex == path.get(i + 1)){
                        v.itsBoundary.add(new Border(Integer.toString(j), "S", "N"));
                    }
                    else if (v.southNeighbor != null && v.southNeighbor.direction != null
                            && !v.southNeighbor.direction.isEmpty() && v.southNeighbor.toVertex == path.get(i + 1)){
                        v.itsBoundary.add(new Border(Integer.toString(j), "E", "W"));
                    }
                    else if (v.eastNeighbor != null && v.eastNeighbor.direction != null
                            && !v.eastNeighbor.direction.isEmpty() && v.eastNeighbor.toVertex == path.get(i + 1)){
                        v.itsBoundary.add(new Border(Integer.toString(j), "N", "S"));
                    }
                    i++;
                }
            }
            j++;
            i = 0;
        }
    }
    private static <T> boolean listEqualsIgnoreOrder(ArrayList<T> list1, ArrayList<T> list2) {
        return new HashSet<>(list1).equals(new HashSet<>(list2));
    }
    private int markLeafs() {
        int minimum_edge_amount = 5;
        List<Vertex> markedOnes = new ArrayList<>();
        for (Vertex v : vertices) {
            if (v.edgeAmount == 1 && v.mark == 0) {
                v.mark = 1;
                markedOnes.add(v);
                for (Edge neighbor : v.neighbors) {
                    if (neighbor.toVertex.mark == 0) {
                        neighbor.toVertex.edgeAmount -= 1;
                    }
                }
            }
        }
        if (markedOnes.size() == vertices.size() || vertices.size() == 1){
            return -1;
        }
        for (Vertex v : vertices) {
            if (v.edgeAmount <= minimum_edge_amount && v.mark != 1) {
                minimum_edge_amount = v.edgeAmount;
            }
        }
        return minimum_edge_amount;
    }
    private void removePathMarkings() {
        for (Vertex v : vertices) {
            v.pathMark = 0;
        }
    }
    private ArrayList<Vertex> getFaces(Vertex v, Edge start) {
        int get_out = 0;
        if (v.northNeighbor != null && v.northNeighbor.direction != null
                && !v.northNeighbor.direction.isEmpty()
                && v.westNeighbor != null && v.westNeighbor.direction != null
                && !v.westNeighbor.direction.isEmpty()){
        }
        else{
            get_out = 1;
        }
        Vertex transv;
        Vertex notTransv;
        int round = 0;
        int trial = 0;
        int flag = 0;
        List<Vertex> goBack = new ArrayList<Vertex>();
        List<Integer> goBackTrial = new ArrayList<Integer>();

        ArrayList<Vertex> finalPath = new ArrayList<Vertex>();
        finalPath.add(v);
        transv = start.toVertex;
        finalPath.add(transv);
        transv.pathMark = 1;
        int i = 0;
        do {
            if (get_out == 1){
                break;
            }
            notTransv = transv;
            if (trial == 0) {
                if (transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()

                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.mark != 1

                        && transv.eastNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    if (round > 0){
                        goBack.add(transv.eastNeighbor.toVertex);
                        goBackTrial.add(3);
                    }
                    goBack.add(transv.northNeighbor.toVertex);
                    goBackTrial.add(0);
                    transv = transv.westNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 1;
                    transv.pathMark = 1;
                }
                else if (transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.northNeighbor.toVertex);
                    goBackTrial.add(0);
                    transv = transv.westNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 1;
                    transv.pathMark = 1;
                }
                else if (transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    if (round > 0){
                        goBack.add(transv.eastNeighbor.toVertex);
                        goBackTrial.add(3);
                    }
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    if (round > 0){
                        goBack.add(transv.eastNeighbor.toVertex);
                        goBackTrial.add(3);
                    }
                    transv = transv.westNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 1;
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    transv = transv.westNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    trial = 1;
                }
                else if (transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (round > 0 && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1) {
                    transv = transv.eastNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    trial = 3;
                }
            } else if (trial == 1) {
                round++;
                if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()

                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.mark != 1

                        && transv.southNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.westNeighbor.toVertex);
                    goBack.add(transv.northNeighbor.toVertex);
                    goBackTrial.add(1);
                    goBackTrial.add(0);
                    transv = transv.southNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 2;
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()

                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.mark != 1

                        && transv.southNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.westNeighbor.toVertex);
                    goBackTrial.add(1);
                    transv = transv.southNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 2;
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.northNeighbor.toVertex);
                    transv = transv.westNeighbor.toVertex;
                    goBackTrial.add(0);
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()

                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.mark != 1

                        && transv.northNeighbor.toVertex.pathMark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.northNeighbor.toVertex);
                    transv = transv.southNeighbor.toVertex;
                    goBackTrial.add(0);
                    finalPath.add(transv);
                    trial = 2;
                    transv.pathMark = 1;
                }
                else if (transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    transv = transv.southNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    trial = 2;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    transv = transv.westNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    trial = 0;
                }
            } else if (trial == 2) {
                round++;
                if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()

                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.mark != 1

                        && transv.eastNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.westNeighbor.toVertex);
                    goBack.add(transv.southNeighbor.toVertex);
                    goBackTrial.add(1);
                    goBackTrial.add(2);
                    transv = transv.eastNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 3;
                    transv.pathMark = 1;
                }
                else if (transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()

                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.mark != 1

                        && transv.eastNeighbor.toVertex.pathMark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.southNeighbor.toVertex);
                    transv = transv.eastNeighbor.toVertex;
                    goBackTrial.add(2);
                    finalPath.add(transv);
                    trial = 3;
                    transv.pathMark = 1;
                }
                else if (transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()

                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.mark != 1

                        && transv.northNeighbor.toVertex.pathMark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.northNeighbor.toVertex);
                    goBackTrial.add(0);
                    transv = transv.southNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()

                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.mark != 1

                        && transv.eastNeighbor.toVertex.pathMark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.westNeighbor.toVertex);
                    transv = transv.eastNeighbor.toVertex;
                    goBackTrial.add(1);
                    finalPath.add(transv);
                    trial = 3;
                    transv.pathMark = 1;
                }
                else if (transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1) {
                    transv = transv.eastNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 3;
                    transv.pathMark = 1;
                }
                else if (transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    transv = transv.southNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.westNeighbor.toVertex.mark != 1
                        && transv.westNeighbor.toVertex.pathMark != 1) {
                    transv = transv.westNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    trial = 1;
                }
            } else if (trial == 3) {
                round++;
                if (transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()

                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.mark != 1

                        && transv.northNeighbor.toVertex.pathMark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.eastNeighbor.toVertex);
                    goBack.add(transv.southNeighbor.toVertex);
                    goBackTrial.add(3);
                    goBackTrial.add(2);
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 0;
                    transv.pathMark = 1;
                }
                else if (transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()

                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.mark != 1

                        && transv.northNeighbor.toVertex.pathMark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.eastNeighbor.toVertex);
                    goBackTrial.add(3);
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 0;
                    transv.pathMark = 1;
                }
                else if (transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.southNeighbor.toVertex);
                    goBackTrial.add(2);
                    transv = transv.eastNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.westNeighbor != null && transv.westNeighbor.direction != null
                        && !transv.westNeighbor.direction.isEmpty()
                        && transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()

                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.mark != 1

                        && transv.northNeighbor.toVertex.pathMark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    goBack.add(transv.southNeighbor.toVertex);
                    goBackTrial.add(2);
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 0;
                    transv.pathMark = 1;
                }
                else if (transv.northNeighbor != null && transv.northNeighbor.direction != null
                        && !transv.northNeighbor.direction.isEmpty()
                        && transv.northNeighbor.toVertex.mark != 1
                        && transv.northNeighbor.toVertex.pathMark != 1) {
                    transv = transv.northNeighbor.toVertex;
                    finalPath.add(transv);
                    trial = 0;
                }
                else if (transv.eastNeighbor != null && transv.eastNeighbor.direction != null
                        && !transv.eastNeighbor.direction.isEmpty()
                        && transv.eastNeighbor.toVertex.mark != 1
                        && transv.eastNeighbor.toVertex.pathMark != 1) {
                    transv = transv.eastNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                }
                else if (transv.southNeighbor != null && transv.southNeighbor.direction != null
                        && !transv.southNeighbor.direction.isEmpty()
                        && transv.southNeighbor.toVertex.mark != 1
                        && transv.southNeighbor.toVertex.pathMark != 1) {
                    transv = transv.southNeighbor.toVertex;
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    trial = 2;
                }
            }
            if (notTransv.equals(transv)){
                if (!goBack.isEmpty()){
                    trial = goBackTrial.get(goBackTrial.size() - 1);
                    transv = goBack.get(goBack.size() - 1);
                    faceFinderReturn(finalPath, transv);
                    finalPath.add(transv);
                    transv.pathMark = 1;
                    goBack.remove(goBack.size() - 1);
                    goBackTrial.remove(goBackTrial.size() - 1);
                    i--;
                }
            }
            i++;
            if (i > vertices.size() * vertices.size()) {
                flag = 1;
                removePathMarkings();
                break;
            }
        } while (!transv.equals(v));
        removePathMarkings();
        if (flag == 1 || get_out == 1) {
            finalPath.clear();
            return finalPath;
        } else {
            return finalPath;
        }
    }
    private void faceFinderReturn(ArrayList<Vertex> finalPath, Vertex transv){
        int i = 0;
        if (!finalPath.isEmpty()) {
            ListIterator<Vertex> ok = finalPath.listIterator(finalPath.size());
            while (ok.hasPrevious()) {
                Vertex holder = ok.previous();
                if (holder.eastNeighbor != null && holder.eastNeighbor.direction != null
                        && !holder.eastNeighbor.direction.isEmpty()
                        && holder.eastNeighbor.toVertex == transv) {
                    break;
                } else if (holder.northNeighbor != null && holder.northNeighbor.direction != null
                        && !holder.northNeighbor.direction.isEmpty()
                        && holder.northNeighbor.toVertex == transv) {
                    break;
                } else if (holder.westNeighbor != null && holder.westNeighbor.direction != null
                        && !holder.westNeighbor.direction.isEmpty()
                        && holder.westNeighbor.toVertex == transv) {
                    break;
                } else if (holder.southNeighbor != null && holder.southNeighbor.direction != null
                        && !holder.southNeighbor.direction.isEmpty()
                        && holder.southNeighbor.toVertex == transv) {
                    break;
                }
                else{
                    i++;
                }
            }
            while (i > 0) {
                finalPath.remove(finalPath.size() - 1);
                i--;
            }
        }
    }
}