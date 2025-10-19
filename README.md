# KruskalWallis

Kruskal-Wallis test with Dunn's post-hoc test and Holm correction for multiple comparisons.

* Kruskal-Wallis Test - A non-parametric test to determine if there are statistically significant differences between groups
* Dunn's Post-Hoc Test - Pairwise comparisons to identify which specific groups differ
* Holm Correction - Adjusts p-values to control for multiple comparisons

## Installation

```elixir
def deps do
  [
    {:kruskal_wallis, "~> 0.1.0"}
  ]
end
```

## Usage

```elixir
groups = %{
  "A" => [1.2, 2.3, 1.8, 2.1, 1.9],
  "B" => [3.4, 4.1, 3.8, 4.5, 3.9],
  "C" => [2.1, 2.8, 2.4, 2.6, 2.3]
}

groups
|> KruskalWallis.test()
|> KruskalWallis.dunn_test(groups)
```
