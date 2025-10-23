defmodule KruskalWallisTest do
  use ExUnit.Case
  doctest KruskalWallis

  test "deterministic ranks regardless of map insertion order" do
    groups1 = %{"A" => [1, 3], "B" => [2, 4]}
    groups2 = %{"B" => [2, 4], "A" => [1, 3]}

    res1 = KruskalWallis.test(groups1)
    res2 = KruskalWallis.test(groups2)

    assert res1.h_statistic == res2.h_statistic
    assert res1.group_ranks == res2.group_ranks
  end

  test "ties produce average ranks" do
    groups = %{"A" => [1, 1], "B" => [1, 2]}
    res = KruskalWallis.test(groups)

    # For values [1,1,1,2] the three 1's should receive the average rank 2.0
    assert res.group_ranks["A"] == [2.0, 2.0]
    assert res.group_ranks["B"] == [2.0, 4.0]
  end

  test "post-hoc tests return expected shape and count" do
    groups = %{"A" => [1, 2], "B" => [3, 4], "C" => [5, 6]}
    kw = KruskalWallis.test(groups)

    dunn = KruskalWallis.Dunn.test(kw, groups)
    conover = KruskalWallis.Conover.test(kw, groups)

    # For 3 groups there should be 3 pairwise comparisons
    assert length(dunn) == 3
    assert length(conover) == 3

    for comp <- dunn ++ conover do
      assert Map.has_key?(comp, :group1)
      assert Map.has_key?(comp, :group2)
      assert Map.has_key?(comp, :p_value)
      assert Map.has_key?(comp, :adjusted_p_value)
      assert Map.has_key?(comp, :significant)
    end
  end

  test "empty groups raise or return a sensible error" do
    groups = %{"A" => [], "B" => [1, 2]}

    assert_raise ArgumentError, fn ->
      KruskalWallis.test(groups)
    end
  end

  test "non-numeric input raises" do
    groups = %{"A" => [1, "x"], "B" => [2, 3]}

    assert_raise ArgumentError, fn ->
      KruskalWallis.test(groups)
    end
  end

  test "very small samples still produce post-hoc comps" do
    groups = %{"A" => [1], "B" => [2], "C" => [3]}
    kw = KruskalWallis.test(groups)

    dunn = KruskalWallis.Dunn.test(kw, groups)
    conover = KruskalWallis.Conover.test(kw, groups)

    assert length(dunn) == 3
    assert length(conover) == 3
  end

  test "all ties lead to near-zero H statistic" do
    groups = %{"A" => [1, 1], "B" => [1, 1], "C" => [1, 1]}
    res = KruskalWallis.test(groups)

    # All values identical => no differences between groups; H should be very small (close to 0)
    assert res.h_statistic < 1.0e-6
  end
end
