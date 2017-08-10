#ifndef GRAPH_H
#define GRAPH_H

typedef struct
{
	int size;
	int max_size;
	void **arr;
	int (*cmp)(void *, void *);
} heap;

typedef struct _sllist
{
	int key;
	struct _sllist *next;
} sllist;

typedef struct
{
	int V, E;
	char **adj;
	double **weight;
} graph;

typedef struct
{
	int u, v;
	double w;
} edge;

typedef struct
{
	int u;
	double w;
} vertex;


/* Vertex utilities */
vertex *create_vertex(int u, double w);
void destroy_vertex(vertex *v);
/* ---------------- */
/* Binary Heap functions */
heap *create_heap(int N, int (*cmp)(void *, void *));
void destroy_heap(heap *h);
int heap_is_empty(heap *h);
void min_heapify(heap *h);
void *heap_min(heap *h);
void *heap_extract_min(heap *h);
int heap_insert(heap *h, void *key);
/* --------------------- */

/* Linked list functions */
sllist *create_sll();
void sll_insert_first(sllist **l, int a);
void sll_insert_last(sllist **l, int a);
sllist *sll_remove_first(sllist **l);
sllist *sll_remove_last(sllist **l);
/* --------------------- */

/* Disjoint Set */
int *create_dj_set(int N);
int dj_set(int set[], int s1);
void dj_union(int set[], int s1, int s2);
/*--------------*/

/* Graph functions */
graph *create_graph(int V);
graph *create_weighted_graph(int V);
int n_vertices(graph *g);
int n_edges(graph *g);
int is_edge(graph *g, int u, int v);
int is_weighted_graph(graph *g);
void destroy_graph(graph *g);
void add_edge(graph *g, int a, int b, double w);
void bfs(graph *g, int s);
void dfs(graph *g, int s);
double kruskal(graph *g, graph **out);
double prim(graph *g, graph **out);
double *dijkstra(graph *g, int node);
/* --------------- */


#endif
