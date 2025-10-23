defmodule KruskalWallis.Dunn do
  @moduledoc """
  Dunn's post-hoc test implementation.

  Public API:
    - test/3: run Dunn post-hoc on a Kruskal-Wallis result and original groups
  """
  alias KruskalWallis.Utils

  @doc """
  Performs Dunn's post-hoc test with Holm correction.

  ## Parameters
    - kw_result: Result from kruskal_wallis/1
    - groups: Original groups map

  ## Returns
    List of pairwise comparison results with adjusted p-values
  """
  @spec test(KruskalWallis.Result.t(), map(), number()) :: [map()]
  def test(kw_result, groups, alpha \\ 0.05) do
    # Flatten groups in the same deterministic (sorted key) order used by `test/2`
    group_names = Map.keys(groups) |> Enum.sort()
    all_values = Enum.flat_map(group_names, &Map.get(groups, &1))
    n = length(all_values)

    group_ranks = kw_result.group_ranks
    all_ranks = Map.get(kw_result, :all_ranks) || Map.values(group_ranks) |> List.flatten()

    # Calculate tie correction factor for Dunn's test
    tie_sum = Utils.calculate_tie_sum(all_ranks)

    # Generate all pairwise comparisons
    comparisons =
      for i <- 0..(length(group_names) - 2),
          j <- (i + 1)..(length(group_names) - 1) do
        g1 = Enum.at(group_names, i)
        g2 = Enum.at(group_names, j)

        ranks1 = Map.get(group_ranks, g1)
        ranks2 = Map.get(group_ranks, g2)

        n1 = length(ranks1)
        n2 = length(ranks2)

        r1_mean = Enum.sum(ranks1) / n1
        r2_mean = Enum.sum(ranks2) / n2

        # Dunn's test standard error with proper tie correction
        # SE = sqrt((N(N+1)/12 - T/(12(N-1))) * (1/n1 + 1/n2))
        # where T = Σ(t³ - t) for tied groups
        variance_base = n * (n + 1) / 12.0

        # Apply tie correction if there are ties
        tie_correction =
          if tie_sum > 0 do
            tie_sum / (12.0 * (n - 1))
          else
            0
          end

        se = :math.sqrt((variance_base - tie_correction) * (1.0 / n1 + 1.0 / n2))

        # Z-statistic
        z = (r1_mean - r2_mean) / se

        # Two-tailed p-value from standard normal
        p_value = 2 * (1 - Utils.standard_normal_cdf(abs(z)))

        %{
          group1: g1,
          group2: g2,
          z_statistic: z,
          p_value: p_value,
          mean_rank1: r1_mean,
          mean_rank2: r2_mean
        }
      end

    Utils.holm_correction(comparisons, alpha)
  end
end
