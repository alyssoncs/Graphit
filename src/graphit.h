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
	size_t V, E;
	char **adj;
	double **weight;
} graph;

typedef struct
{
	size_t u, v;
	double w;
} edge;

typedef struct
{
	size_t u;
	double w;
} vertex;


/* Vertex utilities */
vertex *create_vertex(size_t u, double w);
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
graph *create_graph(size_t V);
graph *create_weighted_graph(size_t V);
size_t n_vertices(graph *g);
size_t n_edges(graph *g);
/* lacks a function to know if there is a edge between two vertices */
int is_weighted_graph(graph *g);
void destroy_graph(graph *g);
void add_edge(graph *g, size_t a, size_t b, double w);
void bfs(graph *g, size_t s);
void dfs(graph *g, size_t s);
double kruskal(graph *g, graph **out);
double prim(graph *g, graph **out);
double *dijkstra(graph *g, size_t node);
/* --------------- */


#endif
