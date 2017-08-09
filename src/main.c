#include <stdio.h>
#include <stdlib.h>
#include "graphit.h"

void print_graph(graph *g)
{
	if (g)
	{
		//printf("n_vertices(g) = %zu\n", n_vertices(g));
		for (int i = 0; i < n_vertices(g); i++)
		{
			for (int j = i+1; j < n_vertices(g); j++)
			{
				//printf("test (%d, %d)\n", i, j);
				if (g->adj[i][j])
					printf("(%d, %d)\n", i, j);
			}
		}
	}
}

int main(int argc, char const *argv[]) {
	size_t v, e, a, b;
	double w;
	graph *out1 = NULL, *out2 = NULL;

	scanf("%zu\n", &v);
	scanf("%zu\n", &e);
	graph *g = create_weighted_graph(v);
	if (g)
	{
		while (e--)
		{
			scanf("%zu %zu %lf\n", &a, &b, &w);
			add_edge(g, a, b, w);
			add_edge(g, b, a, w);
		}


		double w1 = kruskal(g, &out1);
		double w2 = prim(g, &out2);
		//double w1 = kruskal(g, NULL);
		//double w2 = prim(g, NULL);
		double *cost = dijkstra(g, 0);

		printf("\n\tkruskal\n");
		printf("%lf\n", w1);
		if (out1)
		{
			print_graph(out1);
			destroy_graph(out1);
		}

		printf("\n\tprim\n");
		printf("%lf\n", w2);
		if (out2)
		{
			print_graph(out2);
			destroy_graph(out2);
		}

		printf("\n\tdijkstra\n");
		if (cost)
		{
			int i;
			for (i = 0; i < n_vertices(g); i++)
			{
				printf("(0 -> %d) = %lf\n", i, cost[i]);
			}
			printf("\n");
			free(cost);
		}
	}
	destroy_graph(g);
	return 0;
}
