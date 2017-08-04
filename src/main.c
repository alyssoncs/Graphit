#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
/*
int main(int argc, char const *argv[]) {
	heap *h;
	size_t arr[] = {20, 19, 18, 17, 16, 15, 14, 13, 12,
			11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

	h = create_heap(20);
	if (h)
	{
		int i;
		for (i = 0; i < 21; i++)
		{
			if (!heap_insert(h, arr[i]))
			{
				printf("deu merda\n");
			}

		}
		while (!heap_is_empty(h))
		{
			printf("%zu\n", heap_extract_min(h));
		}

	}
	return 0;
}*/

int main(int argc, char const *argv[]) {
	size_t v, s, a, b;
	double w;

	scanf("%zu\n", &v);
	scanf("%zu\n", &s);
	graph *g = create_weighted_graph(v);
	if (g)
	{
		while (scanf("%zu %zu %lf\n", &a, &b, &w) != EOF)
		{
			add_edge(g, a, b, w);
		}

		double *cost = dijkstra(g, s);
		if (cost)
		{
			size_t i;
			for (i = 0; i < g->V; i++)
			{
				printf("{%lf}\n", cost[i]);
			}
			printf("\n");
		}
		else printf("deu ruim\n");
		free(cost);
	}
	destroy_graph(g);
	return 0;
}
