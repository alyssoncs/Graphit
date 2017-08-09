#include <stdio.h>
#include <stdlib.h>
#include "graphit.h"

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
