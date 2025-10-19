# KruskalWallis

Kruskal-Wallis test with post-hoc tests for multiple comparisons.

* Kruskal-Wallis Test - A non-parametric test to determine if there are statistically significant differences between groups
* Dunn's Post-Hoc Test - Pairwise comparisons to identify which specific groups differ. Uses z-statistic (normal distribution) - more conservative, fewer assumptions
* Holm Correction - Adjusts p-values to control for multiple comparisons
* Conover Test - Pairwise comparisons to identify which specific groups differ. Uses t-statistic (t-distribution) - more powerful when distributions have similar shapes, similar to Tukey's HSD for ANOVA. It is generally preferred when you have similar distribution shapes across groups, while Dunn's test is more robust when distributions differ in shape

## Installation

```elixir
def deps do
  [
    {:kruskal_wallis, "~> 0.1.1"}
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

kw_result = KruskalWallis.test(groups)

# Use Dunn test (more conservative)
dunn_results = KruskalWallis.dunn_test(kw_result, groups)

# Or use Conover test (more powerful, same interface!)
conover_results = KruskalWallis.conover_test(kw_result, groups)
```
