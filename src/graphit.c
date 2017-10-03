#include <stdlib.h>
#include <math.h>
#include "graphit.h"

#define heap_parent(i) 	(((i)-1)/2)
#define heap_right(i) 	((2*(i))+2)
#define heap_left(i) 	((2*(i))+1)

/*
 * Struct representing an weighted edge
 *
 * @u 			First vertex of the edge
 * @v 			Second vertex of the edge
 * @w 			The weight of the edge
 */
typedef struct
{
	int u;
	int v;
	double w;
} edge;

/*
 * Struct representing a single vertex with an associated
 * weight, used in Prim's and Dijkstra's algorithms
 * 
 * @u 			The vertex
 * @w 			The weight
 */
typedef struct
{
	int u;
	double w;
} vertex;

/* 
 * Struct for the heap data structure
 *
 * @size 		The current size of the heap
 * @max_size 		The maximum number of elements the heap holds
 * @arr 		The array that will hold the elements
 * @cmp 		A pointer to the  comparison function of the heap.
 * 			The function must return a non-zero if the fist
 * 			argument is less than the second argument
 */
typedef struct
{
	int size;
	int max_size;
	void **arr;
	int (*cmp)(const void *, const void *);
} heap;

/*
 * Struct for the singly linked list data structure
 *
 * @key 		An integer key
 * @next 		A pointer to the next node of the list
 */
typedef struct _sllist
{
	int key;
	struct _sllist *next;
} sllist;


/*
 * Used to sort an array of edges with qsort()
 */
static int cmp_edges(const void *arg1, const void *arg2)
{
	const edge *a = (const edge *)arg1;
	const edge *b = (const edge *)arg2;

	if (a->w > b->w)
		return 1;
	if (a->w < b->w)
		return -1;

	return 0;
}

static vertex *create_vertex(int u, double w)
{
	vertex *v = malloc(sizeof(vertex));

	if (v)
	{
		v->u = u;
		v->w = w;
	}

	return v;
}

static void destroy_vertex(vertex *v)
{
	free(v);
}

/*
 * Used to compare vertices saved in a heap structure.
 */
static int cmp_vertex(const void *arg1, const void *arg2)
{
	const vertex *a = (const vertex *)arg1;
	const vertex *b = (const vertex *)arg2;

	return (a->w < b->w);
}


static void destroy_heap(heap *h)
{
	if (h)
		free(h->arr);

	free(h);
}

/*
* Creates a heap with size N with comparison function cmp
*
* @arg N 		The size of the new heap
* @arg cmp 		The comparison function of the heap
*
* @return 		The address of the new heap, or NULL if it fails
*/
static heap *create_heap(int N, int (*cmp)(const void *, const void *))
{
	heap *h = malloc(sizeof(heap));

	if (h)
	{
		h->arr = malloc(sizeof(void *) * N);
		h->cmp = cmp;
		h->max_size = N;
		h->size = 0;

		if (!h->arr)
		{
			destroy_heap(h);
			return NULL;
		}
	}

	return h;
}

static int is_heap_empty(heap *h)
{
	return (h && h->size == 0);
}

static void _min_heapify_aux(heap *h, int i)
{
	int min;
	int l = heap_left(i);
	int r = heap_right(i);

	if (l < h->size && h->cmp(h->arr[l], h->arr[i]))
		min = l;
	else
		min = i;
	if (r < h->size && h->cmp(h->arr[r], h->arr[min]))
		min = r;

	if (min != i)
	{
		void *tmp       = h->arr[min];
		h->arr[min]     = h->arr[i];
		h->arr[i]       = tmp;

		_min_heapify_aux(h, min);
	}
}

/*
 * Restores min heap property
 *
 * @arg h 		A pointer to the heap to be restored
 */
static void min_heapify(heap *h)
{
	if (h)
		_min_heapify_aux(h, 0);
}

/*
 * Restores min heap property if an element's key was changed
 *
 * @arg h 		A pointer to the heap being updated
 */
static void heap_update(heap *h)
{
	for(int i = h->size/2; i >= 0; i--)
		_min_heapify_aux(h, i);
}

/*
 * Returns the smallest element of a min-heap (the first one)
 * This function will not modify the heap
 *
 * @arg h 		A pointer to the heap
 *
 * @return 		The smallest element of the heap
 */
void *heap_min(heap *h)
{
	if (h)
		return h->arr[0];

	return NULL;
}

/*
 * Removes the smallest element of a min-heap
 *
 * @arg h 		A pointer to the heap
 *
 * @return 		The smallest element of the heap
 */
static void *heap_extract_min(heap *h)
{
	if (h && !is_heap_empty(h))
	{
		void *min = h->arr[0];
		h->arr[0] = h->arr[h->size - 1];
		h->size--;
		min_heapify(h);
		return min;
	}

	return NULL;
}

/*
 * Inserts a new element in the heap
 *
 * @arg h 		A pointer to the heap
 * @arg key 		The new element 
 *
 * @return 		Non-zero in success, zero otherwise
 */
static int heap_insert(heap *h, void *key)
{
	if (h && h->size < h->max_size)
	{
		h->size++;
		h->arr[h->size - 1] = key;

		int i = h->size - 1;
		while (i > 0 && h->cmp(h->arr[i], h->arr[heap_parent(i)]))
		{
			void *tmp               = h->arr[heap_parent(i)];
			h->arr[heap_parent(i)]  = h->arr[i];
			h->arr[i]               = tmp;

			i = heap_parent(i);
		}
		return 1;
	}

	return 0;
}


/*
 * Linked lists were being used in BFS algorithm, but I couldn't find a way to
 * create BFS and DFS functions that would be generic enough to the end user.
 */

/*
static sllist *sll_init(void)
{
	return NULL;
}

static void sll_insert_first(sllist **l, int a)
{
	sllist *node = malloc(sizeof(sllist));

	if (node)
	{
		node->key = a;
		node->next = *l;
		*l = node;
	}
}


static void sll_insert_last(sllist **l, int a)
{
	sllist *node = malloc(sizeof(sllist));

	if (node)
	{
		while (*l)
			l = &(*l)->next;
		node->key = a;
		node->next = *l;
		*l = node;
	}
}

static sllist *sll_remove_first(sllist **l)
{
	sllist *node = *l;

	if (*l)
		*l = node->next;

	return node;
}


static sllist *sll_remove_last(sllist **l)
{
	sllist *node = *l;

	while (*l)
	{
		if (!(*l)->next)
		{
			node = *l;
			*l = (*l)->next;
			break;
		}
		l = &(*l)->next;
	}

	return node;
}
*/

static int *create_dj_set(int N)
{
	int *set = malloc(sizeof(int)*N);

	if (set)
	{
		for (int i = 0; i < N; i++)
			set[i] = i;
	}

	return set;
}

/*
 * Returns the set of element s1 in the disjoint set
 *
 * @param set 		The disjoint set
 * @param s1 		The element that the function will search
 *
 * @return 		An integer representing the s1's set
 */
static int dj_set(int set[], int s1)
{
	if (set[s1] == s1)
		return s1;
	else
		return set[s1] = dj_set(set, set[s1]);
}

/*
 * Performs a union operation of the sets containing s1 and s2
 *
 * @param set 		The disjoint set
 * @param s1 		An element of the first set
 * @param s2 		An element of the second set
 */
static void dj_union(int set[], int s1, int s2)
{
	set[dj_set(set, s1)] = dj_set(set, s2);
}

graph *create_graph(int V)
{
	graph *g = malloc(sizeof(graph));

	if (g)
	{
		g->V = V;
		g->E = 0;
		g->weight = NULL;
		g->adj = malloc(sizeof(char *) * V);
		if (g->adj)
		{
			for (int i = 0; i < V; i++)
				g->adj[i] = NULL;

			for (int i = 0; i < V; i++)
			{
				g->adj[i] = calloc(V, sizeof(char));
				if (!g->adj[i])
				{
					destroy_graph(g);
					g = NULL;
					break;
				}
			}
		}
		else
		{
			destroy_graph(g);
			g = NULL;
		}
	}

	return g;
}

graph *create_weighted_graph(int V)
{
	graph *g = create_graph(V);

	if (g)
	{
		g->weight = malloc(sizeof(double *) * V);
		if (g->weight)
		{
			for (int i = 0; i < V; i++)
				g->weight[i] = NULL;

			for (int i = 0; i < V; i++)
			{
				g->weight[i] = malloc(sizeof(double)*V);
				if (g->weight[i])
				{
					for (int j = 0; j < V; j++)
						g->weight[i][j] = 0.0;
				}
				else
				{
					destroy_graph(g);
					g = NULL;
					break;
				}
			}
		}
		else
		{
			destroy_graph(g);
			g = NULL;
		}
	}

	return g;
}

int is_weighted_graph(graph *g)
{
	return (g && g->weight);
}

/*
 * returns the number of vertices in the graph g
 */
int n_vertices(graph *g)
{
	if (g)
		return g->V;

	return 0;
}

/*
 * returns the number of edges in the graph g
 */
int n_edges(graph *g)
{
	if (g)
		return g->E;

	return 0;
}

/*
 * returns a non-zero if u and v are adjacent vertices,
 * zero otherwise
 */
int is_edge(graph *g, int u, int v)
{
	if (g)
		return g->adj[u][v];

	return 0;
}

/*
 * returns the weight of the edge (u, v)
 */
double edge_weight(graph *g, int u, int v)
{
	if (is_weighted_graph(g))
		return g->weight[u][v];

	return 0.0;
}

/*
 * adds an edge between vertices u and v
 * if g is a weighted graph, adds a weight w to this new edge
 */
void add_edge(graph *g, int u, int v, double w)
{
	if (g)
	{
		g->adj[u][v] = 1;
		g->E++;
		if (is_weighted_graph(g))
			g->weight[u][v] = w;
	}
}

void destroy_graph(graph *g)
{
	if (g)
	{
		if (g->adj)
		{
			for (int i = 0; i < g->V; i++)
				free(g->adj[i]);
		}
		if (g->weight)
		{
			for (int i = 0; i < g->V; i++)
				free(g->weight[i]);
		}
		free(g->adj);
		free(g->weight);
		free(g);
	}
}


/*
void bfs(graph *g, int s)
{
	if (g && g->V > s)
	{
		int *distance 		= malloc(sizeof(int) * g->V);
		int *visited 		= calloc(g->V, sizeof(int));
		sllist *queue 		= sll_init();

		if (visited && distance)
		{
			visited[s] = 1;
			distance[s] = 0;

			sll_insert_first(&queue, s);

			sllist *node;
			while ((node = sll_remove_last(&queue)))
			{
				int i = node->key;
				for (int j = 0; j < g->V; j++)
				{
					if (is_edge(g, i, j) && !visited[j])
					{
						visited[j] = 1;
						distance[j] = distance[i] + 1;
						sll_insert_first(&queue, j);
					}
				}
				free(node);
			}
		}
		free(visited);
		free(distance);
	}
}

static void _dfs_visit(graph *g, int s, int visited[])
{
	visited[s] = 1;

	for (int i = 0; i < g->V; i++)
		if (is_edge(g, s, i) && !visited[i])
			_dfs_visit(g, i, visited);
}

void dfs(graph *g, int s)
{
	if (g)
	{
		int *visited = calloc(g->V, sizeof(int));
		if (visited)
			_dfs_visit(g, s, visited);

		free(visited);
	}
}
*/


/*
 * returns the weight of the minimum-spanning-tree using kruskal's algorithm,
 * if out is an address of a graph writes the MST on it,
 * if out is NULL, does nothing
 */
double kruskal(graph *g, graph **out)
{
	double sum = 0;

	if (out)
	{
		*out = create_weighted_graph(g->V);
		if (!*out)
			return sum;
	}

	if (is_weighted_graph(g))
	{
		edge *A 	= malloc(sizeof(edge) * g->E);
		int *set 	= create_dj_set(g->E);
		if (A && set)
		{
			int count = 0;
			for (int i = 0; i < g->V; i++)
			{
				for (int j = 0; j < g->V; j++)
				{
					if (is_edge(g, i, j))
					{
						A[count].u = i;
						A[count].v = j;
						A[count].w = edge_weight(g, i, j);
						count++;
					}
				}
			}

			qsort(A, g->E, sizeof(edge), cmp_edges);

			for (int i = 0; i < g->E; i++)
			{
				int u = A[i].u;
				int v = A[i].v;
				double w = A[i].w;
				if (dj_set(set, u) != dj_set(set, v))
				{
					if (out)
					{
						add_edge(*out, u, v, w);
						add_edge(*out, v, u, w);
					}
					dj_union(set, u, v);
					sum += w;
				}
			}


		}
		free(A);
		free(set);
	}

	return sum;
}

/*
 * same as kruskal(), but using prim's algorithm
 */
double prim(graph *g, graph **out)
{
	double sum = 0.0;

	if (is_weighted_graph(g))
	{
		vertex **vertices 	= malloc(sizeof(vertex *) * g->V);
		double *cost		= malloc(sizeof(double) * g->V);
		int *parent		= malloc(sizeof(int) * g->V);
		int *visited 		= calloc(g->V, sizeof(int));
		heap *pq		= create_heap(g->V, cmp_vertex);

		if (visited && pq && cost && parent && vertices)
		{
			for (int i = 1; i < g->V; i++)
				cost[i] = INFINITY;
			cost[0] = 0.0;

			for (int i = 0; i < g->V; i++)
			{
				vertices[i] = create_vertex(i, cost[i]);
				if (!vertices[i])
				{
					while (!is_heap_empty(pq))
					{
						vertex *v;
						v = heap_extract_min(pq);
						destroy_vertex(v);
					}
					/* would it be better to use a goto? */
					free(cost);
					free(parent);
					free(visited);
					free(vertices);
					destroy_heap(pq);
					return sum;
				}
				heap_insert(pq, vertices[i]);
			}
			while (!is_heap_empty(pq))
			{
				vertex *w_vertex;
				w_vertex = (vertex *)heap_extract_min(pq);
				int u = w_vertex->u;
				destroy_vertex(w_vertex);

				visited[u] = 1;

				for (int v = 0; v < g->V; v++)
				{
					if (is_edge(g, u, v) && edge_weight(g, u, v) <
						cost[v] && !visited[v])
					{
						cost[v] = edge_weight(g, u, v);
						vertices[v]->w = cost[v];
						parent[v] = u;

						heap_update(pq);
					}
				}

			}

			for (int i = 1; i < g->V; i++)
				sum += cost[i];
			if (out)
			{
				*out = create_weighted_graph(g->V);
				if (*out)
				{
					int v;
					for (int u = 1; u < g->V; u++)
					{
						v = parent[u];
						add_edge(*out, u, v, cost[v]);
						add_edge(*out, v, u, cost[v]);
					}
				}
			}
		}
		free(cost);
		free(parent);
		free(visited);
		free(vertices);
		destroy_heap(pq);
	}

	return sum;
}

/*
 * find the shortest path between node and all vertices
 * and return those values as an double array with size g->V
 */
double *dijkstra(graph *g, int node)
{
	double *cost = malloc(sizeof(double) * g->V);

	if (is_weighted_graph(g))
	{
		int *visited 		= calloc(g->V, sizeof(int));
		heap *pq 		= create_heap(g->V, cmp_vertex);
		vertex *w_vertex 	= create_vertex(node, 0.0);

		if (cost && visited && pq && w_vertex)
		{
			for (int i = 0; i < g->V; i++)
				cost[i] = INFINITY;
			cost[node] = 0.0;

			heap_insert(pq, (void *)w_vertex);
			while (!is_heap_empty(pq))
			{
				w_vertex = (vertex *)heap_extract_min(pq);
				int u = w_vertex->u;
				destroy_vertex(w_vertex);

				if (visited[u])
					continue;
				visited[u] = 1;

				for (int v = 0; v < g->V; v++)
				{
					double w = edge_weight(g, u, v);
					if (is_edge(g, u, v) && cost[u]+w < cost[v])
					{
						cost[v] = cost[u] + w;
						w_vertex = create_vertex(v, cost[v]);
						if (w_vertex)
							heap_insert(pq, (void *)w_vertex);
						else
						{
							while (!is_heap_empty(pq))
							{
								w_vertex = heap_extract_min(pq);
								destroy_vertex(w_vertex);
							}
							free(cost);
							cost = NULL;
						}
					}
				}
			}
		}
		else
		{
			destroy_vertex(w_vertex);
			free(cost);
			cost = NULL;
		}
		free(visited);
		destroy_heap(pq);
	}

	return cost;
}
