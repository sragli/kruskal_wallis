defmodule KruskalWallis do
  @moduledoc """
  Implements Kruskal-Wallis test with Dunn's post-hoc test and Holm correction.
  """

  @doc """
  Performs Kruskal-Wallis test on grouped data.

  ## Parameters
    - groups: Map of group names to lists of values, e.g., %{"A" => [1,2,3], "B" => [4,5,6]}

  ## Returns
    Map with :h_statistic, :p_value, :df, and :significant keys
  """
  def test(groups, alpha \\ 0.05) do
    all_values = groups |> Map.values() |> List.flatten()
    n = length(all_values)

    ranks = rank_values(all_values)
    group_ranks = assign_ranks_to_groups(groups, ranks)

    # Calculate H statistic
    sum_term =
      Enum.reduce(group_ranks, 0, fn {_group, ranks}, acc ->
        n_i = length(ranks)
        r_i = Enum.sum(ranks) / n_i
        acc + n_i * :math.pow(r_i, 2)
      end)

    h = 12 / (n * (n + 1)) * sum_term - 3 * (n + 1)

    # Degrees of freedom
    df = map_size(groups) - 1

    # Chi-square approximation for p-value
    p_value = 1 - chi_square_cdf(h, df)

    %{
      h_statistic: h,
      p_value: p_value,
      df: df,
      significant: p_value < alpha,
      group_ranks: group_ranks
    }
  end

  @doc """
  Performs Dunn's post-hoc test with Holm correction.

  ## Parameters
    - kw_result: Result from kruskal_wallis/1
    - groups: Original groups map

  ## Returns
    List of pairwise comparison results with adjusted p-values
  """
  def dunn_test(kw_result, groups, alpha \\ 0.05) do
    all_values = groups |> Map.values() |> List.flatten()
    n = length(all_values)

    group_names = Map.keys(groups)
    group_ranks = kw_result.group_ranks

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

        # Dunn's test statistic with tie correction
        # Calculate sum of tied ranks cubed
        all_ranks = Map.values(group_ranks) |> List.flatten()
        tie_correction = calculate_tie_correction(all_ranks)

        # Standard error for Dunn's test
        se =
          :math.sqrt(
            (n * (n + 1) / 12.0 - tie_correction / (12.0 * (n - 1))) * (1.0 / n1 + 1.0 / n2)
          )

        z = (r1_mean - r2_mean) / se

        # Two-tailed p-value from standard normal
        p_value = 2 * (1 - standard_normal_cdf(abs(z)))

        %{
          group1: g1,
          group2: g2,
          z_statistic: z,
          p_value: p_value,
          mean_rank1: r1_mean,
          mean_rank2: r2_mean
        }
      end

    holm_correction(comparisons, alpha)
  end

  # Rank values (average rank for ties)
  defp rank_values(values) do
    indexed = Enum.with_index(values)
    sorted = Enum.sort_by(indexed, fn {v, _} -> v end)

    # Group by value to handle ties
    grouped = Enum.group_by(sorted, fn {v, _} -> v end, fn {_, idx} -> idx end)

    # Assign ranks
    {rank_map, _} =
      Enum.reduce(sorted, {%{}, 1}, fn {val, _}, {map, current_rank} ->
        # Check if we've already assigned this value
        if Map.has_key?(map, Enum.at(grouped[val], 0)) do
          {map, current_rank}
        else
          indices = grouped[val]
          count = length(indices)
          # Average rank for ties
          avg_rank = (current_rank + current_rank + count - 1) / 2

          new_map =
            Enum.reduce(indices, map, fn idx, m ->
              Map.put(m, idx, avg_rank)
            end)

          {new_map, current_rank + count}
        end
      end)

    rank_map
  end

  defp assign_ranks_to_groups(groups, rank_map) do
    # Flatten groups in order to match rank_map indices
    {group_ranks, _} =
      Enum.reduce(groups, {%{}, 0}, fn {group_name, values}, {acc, offset} ->
        ranks =
          Enum.map(0..(length(values) - 1), fn i ->
            Map.get(rank_map, offset + i)
          end)

        {Map.put(acc, group_name, ranks), offset + length(values)}
      end)

    group_ranks
  end

  # Calculate tie correction factor
  defp calculate_tie_correction(ranks) do
    # Group ranks by value to find ties
    rank_groups = Enum.group_by(ranks, & &1)

    Enum.reduce(rank_groups, 0, fn {_rank, group}, acc ->
      t = length(group)

      if t > 1 do
        acc + (t * t * t - t)
      else
        acc
      end
    end)
  end

  # Holm correction for multiple comparisons
  defp holm_correction(comparisons, alpha) do
    sorted = Enum.sort_by(comparisons, & &1.p_value)
    m = length(comparisons)

    Enum.with_index(sorted, 1)
    |> Enum.map(fn {comp, i} ->
      adjusted_p = min(comp.p_value * (m - i + 1), 1.0)

      Map.put(comp, :adjusted_p_value, adjusted_p)
      |> Map.put(:significant, adjusted_p < alpha)
    end)
  end

  # Chi-square CDF using regularized incomplete gamma function
  defp chi_square_cdf(x, _df) when x <= 0, do: 0.0

  defp chi_square_cdf(x, df) do
    # P(X <= x) for chi-square with df degrees of freedom
    # This is the regularized lower incomplete gamma function: P(df/2, x/2)
    igam(df / 2, x / 2)
  end

  # Regularized lower incomplete gamma function P(a, x) = γ(a,x) / Γ(a)
  defp igam(_a, x) when x <= 0, do: 0.0

  defp igam(a, x) when x < a + 1 do
    # Use series representation for x < a + 1
    gamma_series(a, x)
  end

  defp igam(a, x) do
    # Use continued fraction for x >= a + 1
    1.0 - gamma_continued_fraction(a, x)
  end

  # Series expansion: P(a,x) = e^(-x) * x^a * Σ(Γ(a)/Γ(a+1+n) * x^n)
  defp gamma_series(a, x, max_iter \\ 200) do
    log_term = -x + a * safe_log(x) - log_gamma(a)
    # Prevent overflow
    exp_term = if log_term < -700, do: 0.0, else: :math.exp(log_term)

    {sum, _} =
      Enum.reduce_while(1..max_iter, {1.0, 1.0}, fn n, {sum, term} ->
        term = term * x / (a + n)
        sum = sum + term

        if abs(term) < abs(sum) * 1.0e-15 do
          {:halt, {sum, term}}
        else
          {:cont, {sum, term}}
        end
      end)

    exp_term * sum
  end

  # Continued fraction: Q(a,x) = e^(-x) * x^a * (1/(x+1-a+ 1*(1-a)/(x+3-a+ 2*(2-a)/(x+5-a+ ...))))
  defp gamma_continued_fraction(a, x, max_iter \\ 200) do
    log_term = -x + a * safe_log(x) - log_gamma(a)
    # Prevent overflow
    exp_term = if log_term < -700, do: 0.0, else: :math.exp(log_term)

    # Compute continued fraction
    cf_result = continued_fraction_compute(a, x, max_iter)
    exp_term / cf_result
  end

  defp continued_fraction_compute(a, x, max_iter) do
    Enum.reduce(1..max_iter, {1.0, x + 1.0 - a, x + 1.0 - a}, fn n, {c, d, result} ->
      an = -n * (n - a)
      b = x + (2 * n + 1) - a

      d = b + an / d
      d = if abs(d) < 1.0e-30, do: if(d >= 0, do: 1.0e-30, else: -1.0e-30), else: d
      c = b + an / c
      c = if abs(c) < 1.0e-30, do: if(c >= 0, do: 1.0e-30, else: -1.0e-30), else: c

      delta = c / d
      result = result * delta

      {c, d, result}
    end)
    |> elem(2)
  end

  defp safe_log(x) when x <= 0, do: -700.0
  defp safe_log(x), do: :math.log(x)

  # Log gamma function using Lanczos approximation
  defp log_gamma(z) when z < 0.5 do
    # Reflection formula
    :math.log(:math.pi()) - :math.log(:math.sin(:math.pi() * z)) - log_gamma(1.0 - z)
  end

  defp log_gamma(z) do
    # Lanczos approximation
    g = 7

    coef = [
      0.99999999999980993,
      676.5203681218851,
      -1259.1392167224028,
      771.32342877765313,
      -176.61502916214059,
      12.507343278686905,
      -0.13857109526572012,
      9.9843695780195716e-6,
      1.5056327351493116e-7
    ]

    z = z - 1.0
    x = Enum.at(coef, 0)

    x =
      Enum.reduce(1..8, x, fn i, acc ->
        acc + Enum.at(coef, i) / (z + i)
      end)

    t = z + g + 0.5
    0.5 * :math.log(2 * :math.pi()) + (z + 0.5) * :math.log(t) - t + :math.log(x)
  end

  # Standard normal CDF approximation
  defp standard_normal_cdf(x) do
    0.5 * (1 + erf(x / :math.sqrt(2)))
  end

  # Error function approximation
  defp erf(x) do
    sign = if x >= 0, do: 1, else: -1
    x = abs(x)

    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * :math.exp(-x * x)

    sign * y
  end
end
