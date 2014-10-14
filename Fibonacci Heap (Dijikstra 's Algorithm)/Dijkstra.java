import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/*Dijkstra is the main class that contains the function for implementing 
 Dijkstra algorithm, Dijkstra algorithm using simple scheme and Dijkstra algorithm
using Fibonacci heap*/
public class Dijkstra {
    public int source,vertices,edges;
    private static final int MAX = 5200;// Max is variable used to define max capacity for graph.
    

    /**
     * The Edge class represents the adjacency list of a vertex. Each edge
     * has an end point and weight associated with it
     */
    public static class Edge {
        public int endPoint; 
        public int Weight; 

        public Edge(int point, int weight) { // the constructor of edge class
            setEP(point);
            setWeight(weight);
        }

        public void setEP(int e) { // set the end point
            endPoint = e;
        }

        public void setWeight(int weight) { // set the weight of the
            // edge
            Weight = weight;
        }

        public int getEndPoint() { // returns the endpoint of the edge
            return endPoint;
        }

        public int getWeight() { //creturns the weight of the edge
            return Weight;
        }
    }
    
    /*
     implemented a new class ArrayList<A> that emulates the properties of ArrayList
     * and to make user defined function as per my needs in it.
     */
    public static class ArrayList<A> {
        private int size = 0;
        
        private Object obj[];
        
        public ArrayList() {
            obj = new Object[MAX];
        }
        
        public void add(A e) {
            if (size == obj.length) {
                verifysize();
            }
            obj[size++] = e;
        }
        
        public void nulify() {
            
            for (int i = 0; i < size; i++) {
                obj[i] = null;
            }
            size = 0;
        }
        
        public int posOf(Object o) {
            if (o == null) {
                for (int i = 0; i < size; i++)
                    if (obj[i] == null)
                        return i;
            } else {
                for (int i = 0; i < size; i++)
                    if (o.equals(obj[i]))
                        return i;
            }
            return -1;
        }
        
        public int size() {
            return size;
        }
        private void verifysize() {
            int newSize = obj.length * 2;
            obj = Arrays.copyOf(obj, newSize);
        }
        
        public A
        set(int pos, A e) {
            if (pos >= 0 && pos < this.size) {
                @SuppressWarnings("unchecked")
                A oValue = (A) ArrayList.this.obj[pos];
                ArrayList.this.obj[pos] = e;
                return oValue;
            }
            return null;
        }
        
        @SuppressWarnings("unchecked")
        public A get(int i) {
            if (i >= obj.length) {
                throw new IndexOutOfBoundsException("Index: " + i + ", Size: "
                                                    + i);
            }
            return (A) obj[i];
        }
    }


    /**
     * It is the class containing the methods for implementing
     * Dijkstra's algorithm using simple scheme
     */
    public static class SimpleDijkstra {
        public final static int Undefined = -1; // declaring a node as undefined
     
        public static int nodecount = 0; 
        
        public static int start; 
        
        public static boolean flag = false; 
        // to check connectivity of the graph

        /**
         * creating a graph by making a list of adjacency lists (List of
         * Class edge). If a node has been visited in the Dijkstra's algorithm (boolean value 1) or
         * still not visited (boolean value of 0)
         */
        final static ArrayList<ArrayList<Edge>> graph = new ArrayList<ArrayList<Edge>>();
        static ArrayList<Integer> dist = new ArrayList<Integer>();
        static ArrayList<Integer> predecessor = new ArrayList<Integer>();
        static ArrayList<Boolean> visited = new ArrayList<Boolean>();

        // function for implementing Dijkstra's algorithm
        public static void dijkstra() {
            initialhelpvector();
            
            dist.set(start, 0); // the distance of the source node is set to 0
            int currentNode = start;
            while (currentNode != Undefined) { 
                
                visited.set(currentNode, true); 
                // visited list of the current node is set to true
                updateNeighbors(graph.get(currentNode), currentNode);
                
                currentNode = ClosestUnvisited(); 
														 
            }
            for (int i = 0; i < dist.size(); i++) { 
                
                if (dist.get(i) > (nodecount * 1000 + 1)) {
                    flag = true;
                }
            }
        }

        /**
         * The initialize function is used to initialize the value of class
         * members to their default values.
         */
        public static void initialhelpvector() {
            for (int i = 0; i < nodecount; i++) {
                dist.add(Integer.MAX_VALUE);
                predecessor.add(Undefined);
                visited.add(false);
            }
        }

        /**
         * This function takes a collection of edges those represent the neighbouring
         * nodes and the node under consideration as its argument and scans
         * through all the nodes in the adjacency list of the node under
         * consideration.
         * */
         
        public static void updateNeighbors(ArrayList<Edge> neighbors,
                                           int node) {
            for (int i = 0; i < neighbors.size(); i++) {
                int neighbor = neighbors.get(i).getEndPoint();
                int weight = neighbors.get(i).getWeight();
                if (dist.get(neighbor) > dist.get(node) + weight) {
                    dist.set(neighbor, dist.get(node) + weight);
                    predecessor.set(neighbor, node);
                }
            }
        }

        /**
         * This function checks the list "visited" and
         * then returns the index of the node that
         * is currently unvisited
         */
        public static int ClosestUnvisited() {
            int a = Integer.MAX_VALUE;
            int y = Undefined;
            for (int i = 0; i < nodecount; i++) {
                if (!visited.get(i) && dist.get(i) < a) {
                    a = dist.get(i);
                    y = i;
                }
            }
            return y;
        }

        /**
         * This is a very crucial function which cleans up all the data present
         * in the various lists so as to prepare them for the next round of
         * Dijkstra's Algorithm
         */
        public static void cleanup() {
            dist.nulify();
            predecessor.nulify();
            visited.nulify();
        }

    }

    // FIBONACCI HEAP

    /**
     * The node class is the class for defining the node structure for
     * implementing fibonacci heap that is implemented using a doubly circular linked 
     * list which has the following fields in it: parent, child, left, right.
     */
    public static class node {
        node child;
        node parent;
        node left;
        node right;
        int key;
        boolean check;
        int degree;

        public node() {
            left = this;
            right = this;
            key = Integer.MAX_VALUE;
            check = false;
        }
    }

    // This class implements the Fibonacci heap data structure
    public static class fibonacci {
        node min; 
        // level trees
        int numnode; 
       
        
        public fibonacci() { // this is the constructor of the fibonacci heap
            nulify();
        } 
        
        public void nulify() { // method to clear all the variable.
          
            min = null;
            numnode = 0;
        }
        
        public boolean isEmpty() { // returns null is the fibonacci heap is empty
            return min == null;
        }
        
        /**
         * This method takes a node and a key as its argument and inserts the
         * node into the fibonacci heap and sets the key of this node to the
         * value present in the argument of the function.
         */
        
        public void insert(node n, int v) {
            n.key = v;
            if (min != null) {
                n.left = min;
                n.right = min.right;
                min.right = n;
                n.right.left = n;
                if (v < min.key) {
                    min = n;
                }
            } else {
                min = n;
            }
            numnode++;
        }
        
        /**
         * This function takes a node and a key as its argument and decreases
         * the key of this node to the value passed in the argument.
         */
        
        public void decreasekey(node n, int k) {
            n.key = k;
            node par = n.parent;
            if ((par != null) && (k < par.key)) {
                cut(n, par);
                cascadingcut(par);
            }
            if (min != null) {
                if (n.key < min.key) {
                    min = n;
                }
            }
        }
        
        /**
        removes the node from the siblings list whose key is decreased.
         * */
        public void cut(node n, node p) {
            n.right.left = n.left;
            n.left.right = n.right;
            p.degree--;
            if (p.child == n) {
                p.child = n.right;
            }
            if (p.degree == 0) {
                p.child = null;
            }
            n.left = min;
            n.right = min.right;
            min.right = n;
            n.right.left = n;
            n.parent = null;
            n.check = false;
        }
        
        /**
         * This function takes a node and checks if the childcut value is true or not.
         */
        public void cascadingcut(node n) {
            node parn = n.parent;
            if (parn != null) {
                if (n.check) {
                    cut(n, parn);
                    cascadingcut(parn);
                } else {
                    n.check = true;
                }
            }
        }
        
        /*
          the "removemin" operation removes the minimum keyed node from the
          fibonacci tree.After removing the minimum node the children are inserted into the top level list. but they are not simply merged but a method called pairwise merging or consolidation is called the "consolidation method is mentioned in the following parts of the code */
       
        public node removemin() {
            node z = min;
            if (z != null) {
                node x = z.child;
                int numchild = z.degree;
                node temp;
                while (numchild > 0) {
                    temp = x.right;
                    
                    x.right.left = x.left;
                    x.left.right = x.right;
                    
                    x.left = min;
                    x.right = min.right;
                    min.right = x;
                    x.right.left = x;
                    
                    x.parent = null;
                    x = temp;
                    numchild--;
                }
                z.right.left = z.left;
                z.left.right = z.right;
                if (z.right == z) {
                    min = null;
                } else {
                    min = z.right;
                    consolidate();
                }
                numnode--;
            }
            return z;
        }
        
        /*
         * This function scan the entire top level list and combine two nodes
         * who have the same degree.
         */
        public void consolidate() {
            int arraysize = numnode + 1;
            node[] array = new node[arraysize];
            for (int i = 0; i < arraysize; i++) {
                array[i] = null;
            }
            
            int rootcount = 0;
            node x = min;
            if (x != null) {
                rootcount++;
                x = x.right;
                while (x != min) {
                    rootcount++;
                    x = x.right;
                }
            }
            
            /**
             * This part of the function scans all the nodes in the top level
             * doubly circular linked list and pairwise combines those nodes who
             * have the same degree by calling the method combine
             */
            
            while (rootcount > 0) {
                int deg = x.degree;
                node next = x.right;
                while (array[deg] != null) {
                    node y = array[deg];
                    if (x.key > y.key) {
                        node temp = x;
                        x = y;
                        y = temp;
                    }
                    link(y, x);
                    array[deg] = null;
                    deg++;
                }
                array[deg] = x;
                x = next;
                rootcount--;
            }
            
            min = null;
            int i=0;
            while(i< arraysize)
            {  
                if (array[i] != null) {
                    if (min != null) {
                        array[i].left.right = array[i].right;
                        array[i].right.left = array[i].left;
                        
                        array[i].left = min;
                        array[i].right = min.right;
                        min.right = array[i];
                        array[i].right.left = array[i];
                        if (array[i].key < min.key) {
                            min = array[i];
                        }
                    } else {
                        min = array[i];
                    }
                } i++;
            }
        }
        
        /**
         * This method combines two nodes whose arguments are two nodes.
         */
        public void link(node n, node par) {
            n.right.left = n.left;
            n.left.right = n.right;
            n.parent = par;
            if (par.child == null) {
                par.child = n;
                n.right = n;
                n.left = n;
            } else {
                node z = par.child;
                n.left = z;
                n.right = z.right;
                z.right = n;
                n.right.left = n;
            }
            par.degree++;
            n.check = false;
        }
    }
    
    


    // this is the class that implements Dijkstra's algorithm using fibonacci
    // heap
    public static class FibonacciDijkstra {
        public final static int Undefined = -1; 
        
        public static int nodecount = 0;
        
        public static int start = 0; 

        /**
         * creating a graph by making a list of adjacency lists (List of
         * Class edge) We declare a list (distance) for maintaining the shortest
         * distance from source node to every other node.
         */
final static ArrayList<ArrayList<Edge>> graph = new ArrayList<ArrayList<Edge>>();
        static ArrayList<Integer> predecessor = new ArrayList<Integer>();
        static ArrayList<Boolean> visited = new ArrayList<Boolean>();
        static ArrayList<node> arrayofnode = new ArrayList<node>();
        static fibonacci dist = new fibonacci();
        public static void dijkstra() {
            initialhelpvector();
            /**
             * We call the initialize function to initialize the value of class
             * members to their default values. 
             */
            dist.decreasekey(arrayofnode.get(start), 0); 
            
            int currentNode = start; 
            
            while (currentNode != Undefined) {
                visited.set(currentNode, true); 
                
                updateNeighbors(graph.get(currentNode), currentNode);
                
                currentNode = ClosestUnvisited(); // this finds the
                
            }
        }

  

        /**
         * This function takes a collection of edges those represent the neighboring
         * nodes and the node under consideration as its argument and then scans
         * through all the nodes in the adjacency list of the node under
         * consideration.
         */
        public static void updateNeighbors(ArrayList<Edge> n,
                                           int node) {
            for (int i = 0; i < n.size(); i++) {
                int next = n.get(i).getEndPoint();
                int cost = n.get(i).getWeight();
                if (arrayofnode.get(next).key > arrayofnode.get(node).key
                        + cost)
                {
                    dist.decreasekey(arrayofnode.get(next),
                            arrayofnode.get(node).key + cost);
                    predecessor.set(next, node);
                }
            }
        }
        

        /**
         * The initialize function is used to initialize the value of class
         * members to their default values. 
         */
        public static void initialhelpvector() {
            for (int i = 0; i < nodecount; i++) {
                node h = new node();
                arrayofnode.add(h);
                dist.insert(arrayofnode.get(i), Integer.MAX_VALUE);
                predecessor.add(Undefined);
                visited.add(false);
            }
        }
        
       

        /**
         * This is a very crucial function which cleans up all the data present
         * in the various lists so as to prepare them for the next round of
         * Dijkstra's Algorithm
         */
        public static void cleanup() {
            predecessor.nulify();
            visited.nulify();
        }
        public static void check(int x,int c)
        {
            if(x==1000&&c==1000)
            {
                System.out.println("HEAP MEMORY WILL GET FULL 'GRAPH NOT CONNECTED'! Please abort and start the server again ");
            }
        }
        public static int ClosestUnvisited() {
            node x = dist.removemin();
            int y = arrayofnode.posOf(x);
            return y;
        }

    }

    // FIBONACCI HEAP CODE


    public static int rread(int node, int d) {
        Random rand = new Random(); 
        int nodes = node;
        
    
        double density = Math.ceil(d/2);// Here I m evaluating max edges in terms of Density ie  n*(n-1)/2*(d/100). d here is already d*n/100.
        
        int max = 0,i=0;
        SimpleDijkstra.graph.nulify();
        FibonacciDijkstra.graph.nulify();
        while(i<nodes) { 
            SimpleDijkstra.graph.add(new ArrayList<Edge>());
            FibonacciDijkstra.graph.add(new ArrayList<Edge>());
            i++;
        }
        for (i = 0; i < (nodes-1); i++) { // loop running till max edges
            for (int j = 0; j < density; j++) {
                int val = rand.nextInt(1000);
                int u = rand.nextInt(nodes);
                int v = rand.nextInt(nodes);
                int ecost = val;
                if (v > u) {
                    if (v > max) {
                        max = v;
                    }
                } else {
                    if (u > max) {
                        max = u;
                    }
                }
                FibonacciDijkstra.graph.get(u).add(new Edge(v, ecost));
                SimpleDijkstra.graph.get(u).add(new Edge(v, ecost));
                
            }
        }
        max = max + 1;
        
        SimpleDijkstra.nodecount = nodes;
        FibonacciDijkstra.nodecount = nodes;
        return max;
    }

    /**
     * This function reads data about graph from file and uses this data to form
     * the graph and the heaps. 
     */
    public static int readfile(String file) throws IOException {
        FileInputStream fstream = new FileInputStream(file); 
       
        DataInputStream in = new DataInputStream(fstream); 
        
        BufferedReader br = new BufferedReader(new InputStreamReader(in)); 
       
        Dijkstra s=new Dijkstra();
        String NewLine;
        int pos,extra=0;
        int counter=1;
        ArrayList<Integer> data = new ArrayList<Integer>();
        while ((NewLine = br.readLine()) != null) { 
            
            SimpleDijkstra.graph.nulify();
            FibonacciDijkstra.graph.nulify();
            StringTokenizer str = new StringTokenizer(NewLine); 
            
            while (str.hasMoreElements()) {
                String l = str.nextToken(); 
                
                pos = Integer.parseInt(l);
                if (counter==1)//To parse source node
                {   counter++;
                    s.source=pos;
                    System.out.println("\nSource is "+pos+"\n");
                }
                else if(counter==2 || extra ==1 ){
                    if (extra ==1)//To parse edges
                    {   s.edges=pos;
                        System.out.println("No of edges is "+ pos+"\n");
                        extra++;}
                    else{
                        s.vertices=pos;//To parse vertices
                    System.out.println("Total no of vertices is "+pos+"\n");
                    extra++;
                    counter++;
                    }
                    
                                   }
                else{
                data.add(pos);
               
                }
            }
        }
        in.close(); 
        int k = data.size();
        int check=s.vertices*(s.vertices-1)/2;
       //Check cases
        if(s.edges>check)
        {
            System.out.println("NO OF EDGES YOU ARE ENTERING IS GREATER THAN POSSIBLE!  EDIT YOUR INPUT FILE ,SORRY ");
            System.exit(0);
            
        }
       if(s.edges!=k/3)
        {
            System.out.println("NO OF EDGES YOU ARE ENTERING IS NOT EQUAL TO EDGES SPECFIED BY YOU IN INPUT FILE!  EDIT YOUR INPUT FILE ,SORRY ");
            System.exit(0);
        }
      
        int max = 0;
        
       
        try{
            
        for (int j = 0; j < k-1; j = j + 3) {
            SimpleDijkstra.graph.add(new ArrayList<Edge>());
            FibonacciDijkstra.graph.add(new ArrayList<Edge>());
   
        }
        for (int j = 0; j < k- 1; j = j + 3) {//Why k/\?
            int u = data.get(j);
            int v = data.get(j + 1);
            if (v > u) {
                if (v > max) {
                    max = v;
                }
            } else {
                if (u > max) {
                    max = u;
                }
            }
            int ecost = data.get(j + 2);
            // These functions add the edge to the graph and in heaps.
            SimpleDijkstra.graph.get(u).add(new Edge(v, ecost));
            SimpleDijkstra.graph.get(v).add(new Edge(u, ecost));
            FibonacciDijkstra.graph.get(u).add(new Edge(v, ecost));
            FibonacciDijkstra.graph.get(v).add(new Edge(u, ecost));
  
        }
        max = max + 1;
        // Update the number of nodes in all the classes.
        SimpleDijkstra.nodecount = max;
        FibonacciDijkstra.nodecount = max;
    }
        catch (NullPointerException e) {
            System.out.print("DISCONNECTED GRAPH ERROR\n");
            System.exit(0);
        }
        
        finally {
        return s.source;
        }
    
    }

    // MAIN FUNCTION

    public static void main(String[] args) throws IOException {
        int num = args.length;  
        
        long start, stop, time;
        SimpleDijkstra.cleanup();
        FibonacciDijkstra.cleanup();
        
        Dijkstra s=new Dijkstra();
        
        if (args[0].equals("-r")) {
            System.out.println("Number of vertices  " + " Density     "
                    + "Simple scheme(msec)   " + "F-heap scheme(msec)");
           

            
                        double percent = Double.parseDouble(args[2]);
                        int source = Integer.parseInt(args[3]);
                         int z = Integer.parseInt(args[1]);
            
                         System.out.println("Loading.........");
            int density = (int)(z*Math.ceil(percent)/100);
                        int j = rread(z, density);
                        start = System.currentTimeMillis(); 
                        
                        if (SimpleDijkstra.flag == false && j == SimpleDijkstra.nodecount) {
                          
                                SimpleDijkstra.start = source;
                                SimpleDijkstra.dijkstra();
                                SimpleDijkstra.cleanup();
                          
                        }
                        Random rand = new Random();
                        stop = System.currentTimeMillis(); 
                          long time_1=stop-start;
                         //Implemenatation DFS
                        while (SimpleDijkstra.flag && j != SimpleDijkstra.nodecount) {
                            j = rread(z, density);
                            System.out.println("DISCONNECTED GRAPH ERROR"); 
                            
                            start = 0;
                            stop = 0;
                            start = System.currentTimeMillis();
                            
                                SimpleDijkstra.start = source;
                                SimpleDijkstra.dijkstra();
                                SimpleDijkstra.cleanup();
                            
                            stop = System.currentTimeMillis();
                        }
                        time = stop - start;
                        System.out.print("\n" + z + "\t\t\t" +  (percent)
                                + "% \t\t" + time_1 + " \t\t    "); 
                       
                        start = 0; 
                        stop = 0; 
                        FibonacciDijkstra.start = source;
                        long time1 = 0;
                        start = System.currentTimeMillis(); 
                        
            
            
                            FibonacciDijkstra.dijkstra();
                            FibonacciDijkstra.arrayofnode.nulify();
                            FibonacciDijkstra.dist.nulify();
                            FibonacciDijkstra.predecessor.nulify();
                            FibonacciDijkstra.visited.nulify();
            
                        stop = System.currentTimeMillis(); 
                        time1 = stop - start;
                        System.out.print(time1 + " \n");
            
        }

        // simple scheme Dijkstra's algorithm in file read mode
        else if (args[0].equals("-s") && num == 2) {
            String argument = args[1];
            s.source=readfile(argument);
            System.out.print("Vertices  ");
            System.out.print("Distance from source" + " ");
            System.out.print("\n");
            for (int source = 0; source < SimpleDijkstra.nodecount; source++) {
                System.out.print(source + " ");
                SimpleDijkstra.start = source;
                try{
                SimpleDijkstra.dijkstra();
                if (SimpleDijkstra.flag) {
                    System.out.println("DISCONNECTED GRAPH ERROR");
                    break;
                } else {
                    for (int i = 0; i < SimpleDijkstra.dist.size(); i++) {
                        if(i==s.source)
                        System.out.print("\t  " + SimpleDijkstra.dist.get(i));//n
                    }
                }
                System.out.print("\n");
                SimpleDijkstra.cleanup();
                }
                catch (NullPointerException e) {
                    System.out.print("DISCONNECTED GRAPH ERROR");
                    break;
                }
                
            }
        }

        // Dijkstra's algorithm using Fibonacci heap in file read mode
        else if (args[0].equals("-f") && num == 2) {
            String argument = args[1];
            s.source=readfile(argument);
            System.out.print("Vertices  ");
            System.out.print("Distance from source" + " ");
            System.out.print("\n");
            for (int source = 0; source < FibonacciDijkstra.nodecount; source++) {
                System.out.print(source + " ");
                FibonacciDijkstra.start = source;
                try {
                    FibonacciDijkstra.dijkstra();
                    for (int i = 0; i < FibonacciDijkstra.arrayofnode.size(); i++) {
                        if(i==s.source)
                        System.out.print("\t  "+ FibonacciDijkstra.arrayofnode.get(i).key);
                    }
                } catch (NullPointerException e) {
                    System.out.print("DISCONNECTED GRAPH ERROR");
                    break;
                }
                System.out.print("\n");
                FibonacciDijkstra.arrayofnode.nulify();
                FibonacciDijkstra.dist.nulify();
                FibonacciDijkstra.predecessor.nulify();
                FibonacciDijkstra.visited.nulify();
            }
        }


    }
}
