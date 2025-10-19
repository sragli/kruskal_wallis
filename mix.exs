defmodule KruskalWallis.MixProject do
  use Mix.Project

  def project do
    [
      app: :kruskal_wallis,
      version: "0.1.0",
      elixir: "~> 1.17",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      name: "KruskalWallis",
      description: description(),
      package: package(),
      source_url: "https://github.com/sragli/kruskal_wallis",
      docs: docs()
    ]
  end

  defp description() do
    "Kruskal-Wallis test with post-hoc tests for multiple comparisons."
  end

  defp package() do
    [
      files: ~w(lib .formatter.exs mix.exs README.md LICENSE CHANGELOG),
      licenses: ["Apache-2.0"],
      links: %{"GitHub" => "https://github.com/sragli/kruskal_wallis"}
    ]
  end

  defp docs() do
    [
      main: "KruskalWallis",
      extras: ["README.md", "LICENSE", "examples.livemd", "CHANGELOG"]
    ]
  end

  def application do
    [
      extra_applications: [:logger]
    ]
  end

  defp deps do
    [
    ]
  end
end
