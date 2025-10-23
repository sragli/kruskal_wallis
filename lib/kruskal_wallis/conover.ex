defmodule KruskalWallis.Conover do
  @moduledoc """
  Conover-Iman post-hoc test implementation.

  Public API:
    - test/3: run Conover post-hoc on a Kruskal-Wallis result and original groups
  """
  alias KruskalWallis.Utils

  @doc """
  Performs Conover-Iman post-hoc test with Holm correction.

  Conover test is generally more powerful than Dunn test but assumes
  the distributions have similar shapes (like Tukey's HSD for ANOVA).

  ## Parameters
    - kw_result: Result from kruskal_wallis/1
    - groups: Original groups map
    @spec test(map(), map(), number()) :: [map()]

  ## Returns
    List of pairwise comparison results with adjusted p-values
  """
  @spec test(KruskalWallis.Result.t(), map()) :: [map()]
  def test(kw_result, groups, alpha \\ 0.05) do
    # Flatten groups in the same deterministic (sorted key) order used by `test/2`
    all_values =
      Map.keys(groups)
      |> Enum.sort()
      |> Enum.flat_map(&Map.get(groups, &1))

    n = length(all_values)
    k = map_size(groups)

    group_names = Map.keys(groups) |> Enum.sort()
    group_ranks = kw_result.group_ranks
    all_ranks = Map.get(kw_result, :all_ranks) || Map.values(group_ranks) |> List.flatten()

    # Calculate pooled variance estimate S²
    # S² = [1/(N-k)] * [Σ(R_i²) - N((N+1)²/4)]
    sum_squared_ranks = Enum.reduce(all_ranks, 0, fn r, acc -> acc + r * r end)

    s_squared =
      if n - k > 0 do
        (sum_squared_ranks - n * :math.pow(n + 1, 2) / 4.0) / (n - k)
      else
        # Degenerate case: no degrees of freedom for pooled variance; use 0 to avoid NaN/Inf
        0.0
      end

    # Calculate tie correction factor
    tie_sum = Utils.calculate_tie_sum(all_ranks)

    # Apply tie correction to variance if there are ties
    s_squared_corrected =
      if tie_sum > 0 and n - k > 0 do
        s_squared * (n - 1 - kw_result.h_statistic) / (n - k)
      else
        s_squared
      end

    # Degrees of freedom for t-distribution
    df = n - k

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

        # Conover test statistic
        # t = (R̄_i - R̄_j) / sqrt(S² * (1/n_i + 1/n_j))
        se = :math.sqrt(s_squared_corrected * (1.0 / n1 + 1.0 / n2))

        # Guard against zero or near-zero standard error (degenerate cases)
        {t_stat, p_value} =
          if se <= 0.0 or se != se do
            {0.0, 1.0}
          else
            t_stat = (r1_mean - r2_mean) / se
            p = 2 * (1 - Utils.t_cdf(abs(t_stat), df))
            {t_stat, p}
          end

        %{
          group1: g1,
          group2: g2,
          t_statistic: t_stat,
          p_value: p_value,
          mean_rank1: r1_mean,
          mean_rank2: r2_mean,
          df: df
        }
      end

    Utils.holm_correction(comparisons, alpha)
  end
end
