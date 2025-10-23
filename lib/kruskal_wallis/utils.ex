defmodule KruskalWallis.Utils do
  @moduledoc """
  Internal numeric utilities used by Kruskal-Wallis and post-hoc tests.

  Exposes:
    - holm_correction/2
    - calculate_tie_sum/1
    - standard_normal_cdf/1
    - chi_square_cdf/2
    - t_cdf/2
  """

  @doc """
  Holm correction for multiple comparisons with monotonicity enforcement.
  """
  @spec holm_correction([map()], number()) :: [map()]
  def holm_correction(comparisons, alpha) do
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

  @doc """
  Calculates tie sum: Σ(t³ - t) for all tied groups.
  """
  @spec calculate_tie_sum([number()]) :: number()
  def calculate_tie_sum(ranks) do
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

  @doc """
  Standard normal CDF approximation.
  """
  @spec standard_normal_cdf(number()) :: number()
  def standard_normal_cdf(x) do
    0.5 * (1 + erf(x / :math.sqrt(2)))
  end

  @doc """
  Chi-square CDF using regularized incomplete gamma function.
  """
  @spec chi_square_cdf(number(), non_neg_integer()) :: number()
  def chi_square_cdf(x, _df) when x <= 0, do: 0.0

  def chi_square_cdf(x, df) do
    # P(X <= x) for chi-square with df degrees of freedom
    # This is the regularized lower incomplete gamma function: P(df/2, x/2)
    igam(df / 2, x / 2)
  end

  @doc """
  Student's t-distribution CDF.
  """
  @spec t_cdf(number(), pos_integer()) :: number()
  def t_cdf(t, df) when t < 0 do
    1.0 - t_cdf(-t, df)
  end

  def t_cdf(t, df) do
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

  defp safe_log(x) when x <= 0, do: -700.0
  defp safe_log(x), do: :math.log(x)

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
