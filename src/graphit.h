#ifndef GRAPH_H
#define GRAPH_H

/*
 * The structure that defines a directed, weighted or unweighted, graph
 *
 * @V: 			Number of vertices
 * @E: 			Number of edges
 * @adj: 		Adjacency matrix
 * @weight 		Weight matrix, holds the weight of every edge, 
 * 			it will be set to NULL if the graph is unweighted
 * 			if the graph is unweighted
 */
typedef struct
{
	int V;
	int E;
	char **adj;
	double **weight;
} graph;


/*
 * Creates a new unweighted graph with V vertices
 *
 * @param V 		Number of vertices
 *
 * @return 		The new graph
 */
graph *create_graph(int V);

/*
 * Create a new weighted graph with V vertices
 *
 * @param V 		Number of vertices
 * 
 * @return 		The new graph
 */
graph *create_weighted_graph(int V);

/*
 * Tells whether a graph is weighted or not
 *
 * @param g 		Pointer to the graph being analized
 *
 * @return 		Non-zero if g is a weighted graph, zero otherwise
 */
int is_weighted_graph(graph *g);

/*
 * Returns the number of vertices of the graph g
 */
int n_vertices(graph *g);

/*
 * Returns the number of edges of the graph g
 */
int n_edges(graph *g);

/* 
 * Tells whether exists an edge between u and v in a given graph or not
 *
 * @param g 		Pointer to the graph being analized
 * @param u 		Integer representing the vertex u
 * @param v 		Integer representing the vertex v
 *
 * @return 		Non-zero if there is an edge between 
 * 			u and v, zero otherwise
 */ 
int is_edge(graph *g, int u, int v);

/*
 * Returns the weight of the edge (u, v)
 * in the weighted graph g.
 * Will return 0.0 if g is not an weighted graph
 * If the edge (u, v) doesent exists the behavior is undefined
 * 
 * @param g 		A pointer to the weighted graph
 * @param u 		Integer representing the vertex u
 * @param v 		Integer representing the vertex v
 * 
 * @return 		The weight between the edge (u, v)
 */
double edge_weight(graph *g, int u, int v);

/*
 * Adds an edge between the vertices u and v with weight w in the graph g,
 * all vertices are zero-based integers.
 * If g is a unweighted graph, the weight w will be ignored
 *
 * @param g 		A pointer to a weighted or unweighted graph
 * @param u 		An integer representing the vertex u in the (u, v) edge
 * @param v 		An integer representing the vertex v in the (u, v) edge
 * @param w 		An double representing the weight of the (u, v) edge,
 * 			it will be ignored if the g is an unweighted graph
 */
void add_edge(graph *g, int u, int v, double w);

/*
 * Remove the edge between vertices u and v in the graph g,
 * it will work whether the edge exists or not
 *
 * @param g 		A pointer to a weighted or unweighted graph
 * @param u 		An integer representing the vertex u in the (u, v) edge
 * @param v 		An integer representing the vertex v in the (u, v) edge
 */
void remove_edge(graph *g, int u, int v);

/*
 * Destroys an instance of a graph (weighted or unweighted)
 *
 * @param g 		A pointer to the graph beeing destroyed
 */
void destroy_graph(graph *g);

/*
void bfs(graph *g, int s);
void dfs(graph *g, int s);
*/


/*
 * Returns the sum of the weights of the minimum spanning tree
 * of the graph g using kruskal's algorithm, if out is an address of a 
 * pointer to a graph, writes the MST on it, if out is NULL, does nothing
 *
 * @param g 		A pointer to a weighted graph that will have its 
 * 			minimum spanning tree calculated
 * @param out 		An address to a pointer to a graph, the MST will 
 * 			be created as a weighted graph, and the pointer
 * 			will point to it. If out is NULL the MST will not
 * 			be created.
 *
 * @return 		The sum of the weights of the MST
 */
double kruskal(graph *g, graph **out);


/*
 * Returns the sum of the weights of the minimum spanning tree
 * of the graph g using prim's algorithm, if out is an address of a 
 * pointer to a graph, writes the MST on it, if out is NULL, does nothing
 *
 * @param g 		A pointer to a weighted graph that will have its 
 * 			minimum spanning tree calculated
 * @param out 		An address to a pointer to a graph, the MST will 
 * 			be created as a weighted graph, and the pointer
 * 			will point to it. If out is NULL the MST will not
 * 			be created.
 *
 * @return 		The sum of the weights of the MST
 */
double prim(graph *g, graph **out);

/*
 * Finds the shortest path between the vertex node to all other 
 * vertices of the weighted graph g using Dijkstra's algorithm
 *
 * @param g 		A pointer to the weighted graph g.
 * @param node 		The initial vertex wich Dijkstra's algorithm
 * 			will search in the graph.
 *
 * @return 		An array with the size equals the number of vertices 
 * 			in the graph g, the value in the element N of this array
 * 			represents the distance of the shortest path from vertex 
 * 			node to vertex N, and the node'th element will be 0.0
 */
double *dijkstra(graph *g, int node);

#endif
