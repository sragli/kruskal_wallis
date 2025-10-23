defmodule KruskalWallis.Result do
  @moduledoc """
  Result struct returned by `KruskalWallis.test/2`.

  Fields:
    - `:h_statistic` - corrected H statistic (number)
    - `:p_value` - p-value from chi-square approximation (number)
    - `:df` - degrees of freedom (non_neg_integer)
    - `:significant` - boolean indicating p < alpha
    - `:group_ranks` - map of group_name => list of ranks
    - `:all_ranks` - flattened list of all ranks
  """

  @type t :: %__MODULE__{
          h_statistic: number(),
          p_value: number(),
          df: non_neg_integer(),
          significant: boolean(),
          group_ranks: map(),
          all_ranks: [number()]
        }

  defstruct h_statistic: 0.0,
            p_value: 1.0,
            df: 0,
            significant: false,
            group_ranks: %{},
            all_ranks: []
end
