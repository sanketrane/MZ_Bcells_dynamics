---
name: matplotlib-legend-cleanup
description: 'Clean up matplotlib legends by removing default black handles/markers, using explicit legend handles, and preserving label text. Use when plot legends show unwanted black symbols or incorrect auto-picked artists.'
argument-hint: 'Describe the plot type and preferred legend style (text-only, colored handles, or mixed).'
user-invocable: true
disable-model-invocation: false
---

# Matplotlib Legend Cleanup

## When to Use
- A legend shows unwanted black markers/lines.
- Legend entries are auto-generated from the wrong artist.
- You want text-only legend annotations without visible handles.
- You need consistent legend behavior across multiple subplot axes.

## Inputs
- Plotting library: matplotlib (optionally seaborn on top of matplotlib axes).
- Axes to fix: single axis or list of subplot axes.
- Style target:
  - Text-only legend (no marker/line swatches).
  - Explicit colored handles matching plotted series.
  - Hybrid (some entries with handles, some annotation-only).

## Procedure
1. Identify how the current legend is being created.
- If using `ax.legend(["label"])`, matplotlib may bind that label to the first plotted artist (often black).
- If relying on implicit labels from seaborn/matplotlib, verify the chosen handles are correct.

2. Choose a legend strategy.
- Use text-only legend for parameter summaries or equations.
- Use explicit handles when the legend should map to plotted groups.

3. Apply one of these patterns.

Text-only legend (no black handle):
```python
from matplotlib.lines import Line2D

text_only = Line2D([], [], linestyle='none')
ax.legend(
    handles=[text_only],
    labels=[r'$\\log_{10}(P)=6.12$'],
    handlelength=0,
    handletextpad=0,
    frameon=False,
)
```

Explicit colored handles:
```python
from matplotlib.lines import Line2D

handles = [
    Line2D([], [], marker='o', linestyle='none', color='#ab2239', label='Host'),
    Line2D([], [], marker='o', linestyle='none', color='#3b9bb3', label='Donor'),
]
ax.legend(handles=handles, frameon=False)
```

4. For subplot grids, standardize legends axis-by-axis.
- Loop over each axis and apply the chosen strategy.
- Keep font size, frame style, and location consistent.

5. Validate visually.
- Confirm black default handles are gone.
- Confirm text remains readable and mathematically formatted.
- Confirm legend does not overlap key data regions.

## Decision Points
- If the legend is annotation-only, prefer text-only dummy handles.
- If legend maps to data groups, use explicit custom handles.
- If seaborn auto-legend conflicts with custom legend, disable/rewrite legend after plotting.

## Completion Checks
- No unwanted black legend symbols remain.
- Legend labels match intended content and order.
- Styling is consistent across all related figures.
- Figure renders without warnings/errors.

## Example Prompts
- `/matplotlib-legend-cleanup remove black symbol from my legend text`
- `/matplotlib-legend-cleanup convert seaborn auto-legend to explicit Host/Donor colors`
- `/matplotlib-legend-cleanup make all subplot legends text-only with no handles`
