# Neighborhood matrix and spatial weights

## Background

Spatial analysis often involves studying how geographical units (e.g., countries, regions, municipalities, or districts) are related to one another. Unlike conventional statistical analysis, spatial analysis accounts for the possibility that observations in nearby locations may be *spatially dependent*. This phenomenon is commonly referred to as [spatial autocorrelation](autocorrelation-cluster.md#spatial-autocorrelation). To represent spatial relationships between geographical units mathematically, a *neighborhood matrix* is constructed. This matrix identifies which units are considered neighbors and provides the basis for defining a *spatial weights matrix*[1][2].

## Neighborhood matrix

A neighborhood matrix describes the spatial relationships between a set of geographical units. It is typically represented as an $$n \cdot n$$ binary matrix, where $$n$$ is the number of spatial units. Each element $$n_{ij}$$ indicates whether spatial unit $$i$$ is a neighbor of spatial unit $$j$$[2][3]:

$$
n_{ij} =
\begin{cases}
1 & \text{if units } i \text{ and } j \text{ are neighbors},\\
0 & \text{otherwise}.
\end{cases}
$$

The diagonal elements in the matrix are conventionally set to zero, since a spatial unit is not considered its own neighbor. Alternatively, the neighborhood matrix can be represented as a vector by excluding the diagonal elements, including $$n \cdot n - 1$$ rows.

Neighborhood relationships can be defined using different criteria, such as geographical distance or contiguity. In the case of *contiguity*, two spatial units are considered neighbors if their geographical boundaries touch, which may be operationalized in several ways[2][3]:
- Rook contiguity: Two units are neighbors if they share a common boundary segment (an edge).
- Queen contiguity: Two units are neighbors if they share either a common boundary segment or at least one vertex (a corner).

## Spatial weights matrix

A spatial weights matrix $$W$$ formalizes the spatial relationships between geographical units by assigning weights to neighboring units. It is derived from the neighborhood matrix and is used in spatial statistical methods, for example, in the calculation of [spatial autocorrelation statistics](autocorrelation-cluster.md) such as [Moran's I](autocorrelation-cluster.md#morans-i).

For binary contiguity-based weights, the spatial weights matrix is initially defined as[2][3]:

$$
w_{ij} =
\begin{cases}
1 & \text{if units } i \text{ and } j \text{ are neighbors},\\
0 & \text{otherwise}.
\end{cases}
$$

Thus, the spatial weights matrix is equivalent to the binary neighborhood matrix before standardization.

A common approach is to apply row standardization, so that the weights in each row sum to one. The standardized weighs are calculated as[2][3]:

$$
w_{ij}^{*} =
\frac{n_{ij}}{\sum_{j=1}^{n} n_{ij}}
$$

where the denominator represents the number of neighbors of spatial unit $$i$$.

Consequently, each neighboring unit receives an equal weight of $$\frac{1}{K}$$. Non-neighboring units receive a weight of zero. For example, if a region has five neighbors, each neighbor has the weighting $$\frac{1}{5}$$. Row standardization ensures that the weights for each spatial unit sum to one, allowing the spatial lag to be interpreted as the average value among its neighboring units. Units without neighbors (so-called islands) require special treatment, since their row sum is zero.


## References

[1] Cliff AD, Ord JK (1981)
*Spatial Processes: Models & Applications*. 
Pion.

[2] Getis A, Aldstadt J (2004) 
Constructing the spatial weights matrix using a local statistic. 
*Geographical Analysis* 36(2), 90–104. [10.1111/j.1538-4632.2004.tb01127.x](https://doi.org/10.1111/j.1538-4632.2004.tb01127.x)

[3] Zhou X, Lin H (2008)
Spatial Weights Matrix. 
In: Shekhar S, Xiong H (eds) *Encyclopedia of GIS*. Springer: Boston, MA. [10.1007/978-0-387-35973-1_1307](https://doi.org/10.1007/978-0-387-35973-1_1307)