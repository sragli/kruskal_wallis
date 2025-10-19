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
  def test(groups, alpha \\ 0.05)

  def test(groups, _alpha) when not is_map(groups) do
    # Basic input validation: groups must be a map of non-empty lists of numbers
    raise ArgumentError, "groups must be a map of group_name => list_of_numbers"
  end

  def test(groups, alpha) do
    # Flatten groups in a deterministic order (sorted keys) so ranks map to group values reliably
    group_keys = Map.keys(groups) |> Enum.sort()

    # Validate group shapes: each key must map to a non-empty list of numbers
    Enum.each(group_keys, fn k ->
      vals = Map.get(groups, k)

      unless is_list(vals) do
        raise ArgumentError, "group #{inspect(k)} must be a list"
      end

      if length(vals) == 0 do
        raise ArgumentError, "group #{inspect(k)} must not be empty"
      end
    end)

    all_values = Enum.flat_map(group_keys, &Map.get(groups, &1))

    unless Enum.all?(all_values, &is_number/1) do
      raise ArgumentError, "all group values must be numeric"
    end

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

    # Tie correction for H statistic
    all_ranks = Map.values(group_ranks) |> List.flatten()
    tie_sum = calculate_tie_sum(all_ranks)

    # Apply tie correction to H if there are ties
    h_corrected =
      if tie_sum > 0 do
        denom = :math.pow(n, 3) - n
        correction = 1 - tie_sum / denom

        # Guard against degenerate cases where correction <= 0 (e.g., all values tied)
        if correction <= 0.0 do
          0.0
        else
          h / correction
        end
      else
        h
      end

    # Degrees of freedom
    df = map_size(groups) - 1

    # Chi-square approximation for p-value
    p_value = 1 - chi_square_cdf(h_corrected, df)

    %{
      h_statistic: h_corrected,
      p_value: p_value,
      df: df,
      significant: p_value < alpha,
      group_ranks: group_ranks,
      all_ranks: all_ranks
    }
  end

  @doc """
  Performs Conover-Iman post-hoc test with Holm correction.

  Conover test is generally more powerful than Dunn test but assumes
  the distributions have similar shapes (like Tukey's HSD for ANOVA).

  ## Parameters
    - kw_result: Result from kruskal_wallis/1
    - groups: Original groups map

  ## Returns
    List of pairwise comparison results with adjusted p-values
  """
  def conover_test(kw_result, groups, alpha \\ 0.05) do
    # Flatten groups in the same deterministic (sorted key) order used by `test/2`
    all_values =
      Map.keys(groups)
      |> Enum.sort()
      |> Enum.flat_map(&Map.get(groups, &1))

    n = length(all_values)
    k = map_size(groups)

    group_names = Map.keys(groups) |> Enum.sort()
    group_ranks = kw_result.group_ranks
    all_ranks = kw_result[:all_ranks] || Map.values(group_ranks) |> List.flatten()

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
    tie_sum = calculate_tie_sum(all_ranks)

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
            p = 2 * (1 - t_cdf(abs(t_stat), df))
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

    holm_correction(comparisons, alpha)
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
    # Flatten groups in the same deterministic (sorted key) order used by `test/2`
    group_names = Map.keys(groups) |> Enum.sort()
    all_values = Enum.flat_map(group_names, &Map.get(groups, &1))
    n = length(all_values)

    group_ranks = kw_result.group_ranks
    all_ranks = kw_result[:all_ranks] || Map.values(group_ranks) |> List.flatten()

    # Calculate tie correction factor for Dunn's test
    tie_sum = calculate_tie_sum(all_ranks)

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
    # Use a deterministic order for group keys (sorted) to match how we flatten values
    group_keys = Map.keys(groups) |> Enum.sort()

    {group_ranks, _} =
      Enum.reduce(group_keys, {%{}, 0}, fn group_name, {acc, offset} ->
        values = Map.get(groups, group_name)

        ranks =
          Enum.map(0..(length(values) - 1), fn i ->
            Map.get(rank_map, offset + i)
          end)

        {Map.put(acc, group_name, ranks), offset + length(values)}
      end)

    group_ranks
  end

  # Calculate tie sum: Σ(t³ - t) for all tied groups
  defp calculate_tie_sum(ranks) do
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

  # Holm correction for multiple comparisons with monotonicity enforcement
  defp holm_correction(comparisons, alpha) do
    sorted = Enum.sort_by(comparisons, & &1.p_value)
    m = length(comparisons)

    # Calculate adjusted p-values with monotonicity
    {adjusted_comps, _} =
      Enum.reduce(Enum.with_index(sorted, 1), {[], 0.0}, fn {comp, i}, {acc, max_p} ->
        raw_adjusted = comp.p_value * (m - i + 1)
        # Ensure monotonicity: adjusted p-value should not be less than previous
        adjusted_p = min(:erlang.max(raw_adjusted, max_p), 1.0)

        updated_comp =
          comp
          |> Map.put(:adjusted_p_value, adjusted_p)
          |> Map.put(:significant, adjusted_p < alpha)

        {[updated_comp | acc], adjusted_p}
      end)

    Enum.reverse(adjusted_comps)
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

  # Student's t-distribution CDF
  defp t_cdf(t, df) when t < 0 do
    1.0 - t_cdf(-t, df)
  end

  defp t_cdf(t, df) do
    # Use the relationship with beta distribution
    # P(T <= t) = 1 - 0.5 * I_x(df/2, 1/2)
    # where x = df/(df + t²) and I is the regularized incomplete beta function
    x = df / (df + t * t)

    if t >= 0 do
      1.0 - 0.5 * incomplete_beta(x, df / 2.0, 0.5)
    else
      0.5 * incomplete_beta(x, df / 2.0, 0.5)
    end
  end

  # Regularized incomplete beta function I_x(a,b)
  defp incomplete_beta(x, _a, _b) when x <= 0, do: 0.0
  defp incomplete_beta(x, _a, _b) when x >= 1, do: 1.0

  defp incomplete_beta(x, a, b) do
    # Use continued fraction representation
    bt =
      :math.exp(
        log_gamma(a + b) - log_gamma(a) - log_gamma(b) +
          a * :math.log(x) + b * :math.log(1.0 - x)
      )

    # Use symmetry if x > (a+1)/(a+b+2)
    if x < (a + 1.0) / (a + b + 2.0) do
      bt * beta_continued_fraction(x, a, b) / a
    else
      1.0 - bt * beta_continued_fraction(1.0 - x, b, a) / b
    end
  end

  # Continued fraction for incomplete beta function
  defp beta_continued_fraction(x, a, b, max_iter \\ 200) do
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    # First iteration
    c = 1.0
    d = 1.0 - qab * x / qap
    d = if abs(d) < 1.0e-30, do: 1.0e-30, else: d
    d = 1.0 / d
    result = d

    {final_result, _} =
      Enum.reduce_while(1..max_iter, {result, {c, d}}, fn m, {h, {c, d}} ->
        m_float = m * 1.0
        m2 = 2.0 * m_float

        # Even iteration
        aa = m_float * (b - m_float) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = if abs(d) < 1.0e-30, do: 1.0e-30, else: d
        c = 1.0 + aa / c
        c = if abs(c) < 1.0e-30, do: 1.0e-30, else: c
        d = 1.0 / d
        h = h * d * c

        # Odd iteration
        aa = -(a + m_float) * (qab + m_float) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = if abs(d) < 1.0e-30, do: 1.0e-30, else: d
        c = 1.0 + aa / c
        c = if abs(c) < 1.0e-30, do: 1.0e-30, else: c
        d = 1.0 / d
        delta = d * c
        h = h * delta

        if abs(delta - 1.0) < 1.0e-15 do
          {:halt, {h, {c, d}}}
        else
          {:cont, {h, {c, d}}}
        end
      end)

    final_result
  end
end
