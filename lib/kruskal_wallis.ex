defmodule KruskalWallis do
  @moduledoc """
  Implements the Kruskal-Wallis test.

  Public API:
    - test/2: Perform Kruskal-Wallis on grouped numeric data.
  """
  alias KruskalWallis.Utils

  @doc """
  Performs Kruskal-Wallis test on grouped data.

  ## Parameters
    - groups: Map of group names to lists of values, e.g., %{"A" => [1,2,3], "B" => [4,5,6]}

  ## Returns
    Map with :h_statistic, :p_value, :df, and :significant keys
  """
  @spec test(map(), number()) :: KruskalWallis.Result.t()
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
    tie_sum = Utils.calculate_tie_sum(all_ranks)

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
    p_value = 1 - Utils.chi_square_cdf(h_corrected, df)

    %KruskalWallis.Result{
      h_statistic: h_corrected,
      p_value: p_value,
      df: df,
      significant: p_value < alpha,
      group_ranks: group_ranks,
      all_ranks: all_ranks
    }
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
end
