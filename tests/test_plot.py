import dionysus as d
from matplotlib.figure import Figure


def test_plot_bars_marks_infinite_interval_beyond_finite_deaths():
    diagram = d.Diagram([
        (0, float("inf")),
        (0, 3),
        (0, 5),
        (0, 3),
    ])
    ax = Figure().subplots()

    d.plot.plot_bars(diagram, ax=ax)

    infinite_lines = [line for line in ax.lines if line.get_marker() == ">"]
    finite_lines = [line for line in ax.lines if line.get_marker() != ">"]
    assert len(infinite_lines) == 1
    assert infinite_lines[0].get_xdata()[-1] > max(
        line.get_xdata()[-1] for line in finite_lines
    )
    assert infinite_lines[0].get_markevery() == [1]


def test_plot_bars_infinite_interval_advances_at_large_magnitude():
    diagram = d.Diagram([(-1e16, float("inf"))])
    ax = Figure().subplots()

    d.plot.plot_bars(diagram, ax=ax)

    xdata = ax.lines[0].get_xdata()
    assert xdata[-1] > xdata[0]
