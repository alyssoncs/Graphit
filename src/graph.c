#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "graph.h"

#define heap_parent(i) 	(size_t)( ((i)-1)/2 )
#define heap_right(i) 	(size_t)( (2*(i))+2 )
#define heap_left(i) 	(size_t)( (2*(i))+1 )


vertex *create_vertex(size_t u, double w)
{
	vertex *v = (vertex *)malloc(sizeof(vertex));
	if (v)
	{
		v->u = u;
		v->w = w;
		return v;
	}
	else return NULL;
}
void destroy_vertex(vertex *v)
{
	free(v);
}


/*
* Creates a heap with size N
* and with comparison function cmp
*/
heap *create_heap(size_t N, int (*cmp)(void *, void *))
{
	heap *h = (heap*)malloc(sizeof(heap));
	if (h)
	{
		h->arr = (void *)malloc(sizeof(void *)*N);
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

void destroy_heap(heap *h)
{
	if (h)
		free(h->arr);

	free(h);
}

int heap_is_empty(heap *h)
{
	if (h)
	{
		if(h->size == 0)
		return 1;
	}
	return 0;
}

static void _min_heapify_aux(heap *h, size_t i)
{
	size_t min;
	size_t l = heap_left(i);
	size_t r = heap_right(i);
	if (l < h->size && h->cmp(h->arr[l], h->arr[i]))
		min = l;
	else
		min = i;
	if (r < h->size && h->cmp(h->arr[r], h->arr[min]))
		min = r;

	if (min != i)
	{
		void *tmp 	= h->arr[min];
		h->arr[min] 	= h->arr[i];
		h->arr[i] 	= tmp;

		_min_heapify_aux(h, min);
	}
}

void min_heapify(heap *h)
{
	if (h)
		_min_heapify_aux(h, 0);
}

void *heap_min(heap *h)
{
	if (h)
		return h->arr[0];
	return NULL;
}

void *heap_extract_min(heap *h)
{
	if (h && !heap_is_empty(h))
	{
		void *min = h->arr[0];
		h->arr[0] = h->arr[h->size-1];
		h->size--;
		min_heapify(h);
		return min;
	}
	return NULL;
}

int heap_insert(heap *h, void *key)
{
	if(h && h->size < h->max_size)
	{
		h->size++;
		h->arr[h->size-1] = key;

		size_t i = h->size-1;
		while (i > 0 && h->cmp(h->arr[i], h->arr[heap_parent(i)]))
		{
			void *tmp 		= h->arr[heap_parent(i)];
			h->arr[heap_parent(i)] 	= h->arr[i];
			h->arr[i] 		= tmp;

			i = heap_parent(i);
		}
		return 1;
	}
	return 0;
}

void sll_insert_first(sllist **l, size_t a)
{
	sllist *node = (sllist*)malloc(sizeof(sllist));
	if (node)
	{
		node->key = a;
		node->next = *l;
		*l = node;
	}
}

void sll_insert_last(sllist **l, size_t a)
{
	sllist *node = (sllist*)malloc(sizeof(sllist));
	if (node)
	{
		while (*l)
			l = &(*l)->next;
		node->key = a;
		node->next = *l;
		*l = node;
	}
}

sllist *sll_remove_first(sllist **l)
{
	sllist *node = *l;
	if (*l)
		*l = node->next;
	return node;
}

sllist *sll_remove_last(sllist **l)
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

size_t dj_set(size_t *set, size_t s1)
{
	if (set[s1] == s1)
		return s1;
	else
		return set[s1] = dj_set(set, set[s1]);
}

void dj_union(size_t *set, size_t s1, size_t s2)
{
	set[dj_set(set, s1)] = dj_set(set, s2);
}

graph *create_graph(size_t V)
{
	graph *g = malloc(sizeof(graph));
	if (!g) return NULL;

	g->V = V;
	g->E = 0;
	g->weight = NULL;
	g->adj = (char**)malloc(sizeof(char*)*V);
	if (g->adj)
	{
		size_t i;
		for (i = 0; i < V; i++)
		{
			g->adj[i] = (char*)calloc(V, sizeof(char));
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
	return g;

}

graph *create_weighted_graph(size_t V)
{
	graph *g = create_graph(V);
	if (g)
	{
		g->weight = (double**)malloc(sizeof(double*)*V);
		if (g->weight)
		{
			size_t i;
			for (i = 0; i < V; i++)
			{
				g->weight[i] = (double*)calloc(V, sizeof(double));
				if (!g->weight[i])
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

void destroy_graph(graph *g)
{
	if (g)
	{
		size_t i;
		for (i = 0; i < g->V; i++)
		{
			free(g->adj[i]);
			if (g->weight)
				free(g->weight[i]);
		}
		free(g->adj);
		free(g->weight);
		free(g);
	}
}

void add_edge(graph *g, size_t a, size_t b, double w)
{
	if (g)
	{
		g->adj[a][b] = 1;
		g->E++;
		if (g->weight)
			g->weight[a][b] = w;
	}
}

void bfs(graph *g, size_t s)
{
	if (g && g->V > s)
	{
		size_t *distance 	= (size_t*)malloc(sizeof(size_t)*g->V);
		int *visited		= (int*)calloc(g->V, sizeof(int));
		sllist *queue 		= NULL;

		if (visited && distance)
		{
			size_t i, j;
			visited[s] = 1;
			distance[s] = 0;

			sll_insert_first(&queue, s);

			sllist *node;
			while ((node = sll_remove_last(&queue)))
			{
				i = node->key;
				for (j = 0; j < g->V; j++)
				{
					if (g->adj[i][j] && !visited[j])
					{
						visited[j] = 1;
						distance[j] = distance[i]+1;
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

static void _dfs_visit(graph *g, size_t s, int visited[])
{
	printf("%zu\n", s);
	visited[s] = 1;
	size_t i;
	for (i = 0; i < g->V; i++)
	{
		if (g->adj[s][i] && !visited[i])
		{
			_dfs_visit(g, i, visited);
		}
	}
}

void dfs(graph *g, size_t s)
{
	if (g)
	{
		int *visited = (int*)calloc(g->V, sizeof(int));
		if (visited)
			_dfs_visit(g, s, visited);

		free(visited);
	}
}

static int cmp_edges(const void *arg1, const void *arg2)
{
	edge a = *(edge *)arg1;
	edge b = *(edge *)arg2;

	if (a.w > b.w) return 1;
	if (a.w < b.w) return -1;
	return 0;
}

static int cmp_vertex(void *arg1, void *arg2)
{
	vertex *a = (vertex *)arg1;
	vertex *b = (vertex *)arg2;

	if (a->w < b->w) return 1;
	else return 0;

}

double kruskal(graph *g, graph **out)
{
	size_t i, j;
	double sum = 0;

	if (g && g->weight)
	{
		*out = create_weighted_graph(g->V);
		edge *A = (edge*)malloc(sizeof(edge)*g->E);
		size_t *set = (size_t*)malloc(sizeof(size_t)*g->E);
		if (*out && A && set)
		{
			size_t count = 0;
			for (i = 0; i < g->V; i++)
			{
				for (j = 0; j < g->V; j++)
				{
					if (g->adj[i][j])
					{
						A[count].u = i;
						A[count].v = j;
						A[count].w = g->weight[i][j];
						count++;
					}
				}
			}

			qsort(A, g->E, sizeof(edge), cmp_edges);

			for (i = 0; i < g->E; i++)
				set[i] = i;

			for (i = 0; i < g->E; i++)
			{
				size_t u = A[i].u;
				size_t v = A[i].v;
				double w = A[i].w;
				if (dj_set(set, u) != dj_set(set, v))
				{
					add_edge(*out, u, v, w);
					add_edge(*out, v, u, w);
					dj_union(set, u, v);
					sum += w;
				}
			}


		}
		else
			free(*out);
		free(A);
		free(set);
	}
	return sum;
}



double *dijkstra(graph *g, size_t node)
{
	double *cost = (double*)malloc(sizeof(double)*g->V);

	if (g && g->weight)
	{

		int *visited = (int*)calloc(g->V, sizeof(int));
		heap *pq = create_heap(g->V, cmp_vertex);
		vertex *w_vertex = create_vertex(node, 0.0);

		if (cost && visited && pq && w_vertex)
		{
			size_t i;
			for (i = 0; i < g->V; i++)
				cost[i] = DBL_MAX;
			cost[0] = 0.0;

			heap_insert(pq, (void *)w_vertex);
			while (!heap_is_empty(pq))
			{
				w_vertex = (vertex *)heap_extract_min(pq);
				size_t u = w_vertex->u;
				free(w_vertex);

				if (visited[u])
					continue;
				visited[u] = 1;

				size_t v;
				for (v = 0; v < g->V; v++)
				{
					if (g->adj[u][v])
					{
						double w = g->weight[u][v];
						if (cost[u]+w < cost[v])
						{
							cost[v] = cost[u]+w;
							w_vertex = create_vertex(v, cost[v]);
							if (w_vertex)
							heap_insert(pq, (void *)w_vertex);
							else return NULL;
						}
					}
				}
			}
			w_vertex = NULL;
		}
		else
		{
			free(visited);
			free(pq);
			free(w_vertex);
			cost = NULL;
		}
	}
	return cost;
}
