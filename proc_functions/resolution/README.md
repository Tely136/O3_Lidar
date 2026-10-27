# Vertical resolution

Leblanc et al. (2016), Part 1, sections 4.1 and 4.2.

```matlab
dzIR = resolution_IR(C,dz);
dzFC = resolution_FC(C,dz);
```

`C` can be the altitude-dependent cell array from `filter_gen`, a single
filter matrix (column 2 is used), or a coefficient vector. Pass column 1
explicitly for an SG smoothing filter. Symmetry identifies smoothing or
first-derivative coefficients, which are normalized internally. `dz` and
the returned resolutions use the same length unit. Cell inputs preserve
their shape. These are filter resolutions; retrieval boundary bins are
not masked.

`resolution_IR` measures the delta response for smoothing or the Heaviside
step response for differentiation, using the outermost half-maximum
crossings and linear interpolation. `resolution_FC` uses the first 0.5
crossing of the signed gain normalized by the ideal filter response.
If no crossing occurs below Nyquist, it returns `dz`.

Optional outputs are `[dzIR,response]` and `[dzFC,gain,f]`, where `f` is in
cycles per bin. For cell inputs, these outputs are also cell arrays.

For successive filters with constant coefficients on the same grid,
convolve the coefficients before computing resolution:

```matlab
c = conv(c_smooth,c_derivative);
dzIR = resolution_IR(c,dz);
dzFC = resolution_FC(c,dz);
```

This supports smoothing followed by at most one first derivative. For
altitude-dependent successive filters, the full composed operator must
be formed using the coefficients at each contributing bin; convolving
only the filters at the central bin is not exact.
