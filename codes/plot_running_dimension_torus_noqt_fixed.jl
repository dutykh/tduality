#!/usr/bin/env julia

# plot_running_dimension_torus_noqt.jl
#
# Self-contained Julia workflow for the potential-based running interaction
# dimension in finite-size T-duality electrostatics on a square torus.
#
# This script avoids Plots.jl, GR.jl, Qt, Makie, and GUI backends. It uses
# SpecialFunctions.jl only for the modified Bessel functions K0 and K1.
# It writes CSV data, a vector PDF figure, and a convergence report.
#
# Run:
#     julia plot_running_dimension_torus_noqt.jl
#
# Formulae implemented for axial separation r=(xi L,0):
#
# U_T2(xi; a) = Delta V_{L,l0}(r,0)/Q^2
#              = sum_{n in Z^2 \ {0}}
#                a*K1(2*pi*a*|n|)/(2*pi*|n|)
#                * [1 - cos(2*pi*n1*xi)],
# where a=l0/L and xi=r/L.
#
# D_T2(xi; a) = 3 - xi*dU_T2/dxi / U_T2.
#
# Infinite-plane reference:
# D_inf(s) = 3 - 2*s^2/((1+s^2)*log(1+s^2)),  s=xi/a.

using Pkg
using Printf
using Dates

const SCRIPT_DIR = @__DIR__
const ENV_DIR = joinpath(SCRIPT_DIR, ".running_dimension_env_noqt")

function ensure_dependencies()
    println("Activating isolated Julia environment: ", ENV_DIR)
    Pkg.activate(ENV_DIR)
    try
        if Base.find_package("SpecialFunctions") === nothing
            Pkg.add("SpecialFunctions")
        end
        Pkg.instantiate()
    catch err
        println(stderr, "Failed to install or instantiate SpecialFunctions.jl.")
        println(stderr, "This script intentionally avoids Plots.jl, GR.jl, Qt, and Makie.")
        rethrow(err)
    end
end

ensure_dependencies()
using SpecialFunctions

# --------------------------- user parameters ---------------------------

const A_VALUES = [0.02, 0.05, 0.10]       # a = l0/L
const XI_MIN = 0.0
const XI_MAX = 0.5
const N_XI = 501

# Potential-tail target used to choose the initial lattice cutoff.
const TAIL_TOL = 1.0e-10

# Direct convergence target for the running dimension curve.
const DIM_CONV_TOL = 1.0e-8

const MAX_N = 4000
const OUT_DIR = joinpath(SCRIPT_DIR, "running_dimension_output_noqt")
const TWOPI = 2.0*pi

# Plot options. The caption in the manuscript can define the y-axis quantity.
const SHOW_TITLE = false
const SHOW_Y_AXIS_LABEL = false

# --------------------------- numerical core ----------------------------

function sanitize_a(a::Float64)::String
    return replace(@sprintf("%.3f", a), "." => "p")
end

# Practical tail estimate for the potential sum. The running dimension is
# checked directly by comparing the N and 2N curves.
function tail_estimate(a::Float64, N::Int)::Float64
    return besselk(0, TWOPI*a*N)/pi
end

function choose_cutoff_from_tail(a::Float64; tol::Float64=TAIL_TOL,
                                 minN::Int=16, maxN::Int=MAX_N)::Int
    N = minN
    while tail_estimate(a, N) > tol
        N *= 2
        if N > maxN
            @warn "Reached maxN before tail target was achieved" a tol N tail=tail_estimate(a, N)
            return maxN
        end
    end

    lo = max(minN, N ÷ 2)
    hi = N
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        if tail_estimate(a, mid) <= tol
            hi = mid
        else
            lo = mid
        end
    end
    return hi
end

# Axial coefficients:
# U(xi;a) = sum_{m=1}^N C_m(a,N) [1 - cos(2*pi*m*xi)].
# The n1=0 modes do not contribute to axial separations, and the +/-m modes
# are combined into the factor 2 below.
function axial_coefficients(a::Float64, N::Int)::Vector{Float64}
    coeffs = zeros(Float64, N)
    for m in 1:N
        acc = 0.0
        for n2 in -N:N
            rho = hypot(m, n2)
            z = TWOPI*a*rho
            acc += a*besselk(1, z)/(TWOPI*rho)
        end
        coeffs[m] = 2.0*acc
    end
    return coeffs
end

function torus_dimension_axial(xis::Vector{Float64}, a::Float64, N::Int)::Vector{Float64}
    coeffs = axial_coefficients(a, N)
    D = zeros(Float64, length(xis))

    for (j, xi) in enumerate(xis)
        if xi == 0.0
            D[j] = 1.0
            continue
        end

        U = 0.0
        xi_dU = 0.0
        for m in 1:N
            qxi = TWOPI*m*xi
            c = coeffs[m]
            U += c*(1.0 - cos(qxi))
            xi_dU += c*qxi*sin(qxi)
        end

        if U <= 1.0e-300
            D[j] = 1.0
        else
            D[j] = 3.0 - xi_dU/U
        end
    end

    return D
end

function infinite_dimension_reference_xi(xi::Float64, a::Float64)::Float64
    s = xi/a
    if s == 0.0
        return 1.0
    elseif abs(s) < 1.0e-5
        # Small-s expansion: 1 + s^2 - 5 s^4/6 + 3 s^6/4 + ...
        s2 = s*s
        return 1.0 + s2 - (5.0/6.0)*s2*s2 + (3.0/4.0)*s2*s2*s2
    else
        return 3.0 - (2.0*s*s)/((1.0 + s*s)*log1p(s*s))
    end
end

function dimension_convergence_error(xis::Vector{Float64}, a::Float64, N::Int)::Float64
    D_N = torus_dimension_axial(xis, a, N)
    D_2N = torus_dimension_axial(xis, a, 2*N)
    return maximum(abs.(D_2N .- D_N))
end

function choose_cutoff_for_dimension(xis::Vector{Float64}, a::Float64)
    N = choose_cutoff_from_tail(a)
    while true
        err = dimension_convergence_error(xis, a, N)
        if err <= DIM_CONV_TOL || 2*N >= MAX_N
            return N, tail_estimate(a, N), err
        end
        N *= 2
    end
end

# -------------------------- PDF plotting code --------------------------

struct Curve
    x::Vector{Float64}
    y::Vector{Float64}
    label::String
    color::Tuple{Float64,Float64,Float64}
    dash::Bool
end

function pdf_escape(s::String)::String
    s = replace(s, "\\" => "\\\\")
    s = replace(s, "(" => "\\(")
    s = replace(s, ")" => "\\)")
    return s
end

function pdf_text(x::Float64, y::Float64, size::Int, text::String; align::Symbol=:left)::String
    # Simple approximate alignment using character count. Helvetica average width ~0.5*font size.
    escaped = pdf_escape(text)
    dx = 0.0
    if align == :center
        dx = -0.25*size*length(text)
    elseif align == :right
        dx = -0.50*size*length(text)
    end
    return @sprintf("BT /F1 %d Tf %.3f %.3f Td (%s) Tj ET\n", size, x+dx, y, escaped)
end

function set_stroke_color(c::Tuple{Float64,Float64,Float64})::String
    return @sprintf("%.3f %.3f %.3f RG\n", c[1], c[2], c[3])
end

function set_fill_color(c::Tuple{Float64,Float64,Float64})::String
    return @sprintf("%.3f %.3f %.3f rg\n", c[1], c[2], c[3])
end

function line_cmd(x1, y1, x2, y2)::String
    return @sprintf("%.3f %.3f m %.3f %.3f l S\n", x1, y1, x2, y2)
end

function polyline_cmd(curve::Curve, xmin, xmax, ymin, ymax, left, bottom, plotW, plotH)::String
    function tx(x)
        return left + (x - xmin)/(xmax - xmin)*plotW
    end
    function ty(y)
        return bottom + (y - ymin)/(ymax - ymin)*plotH
    end

    buf = IOBuffer()
    print(buf, set_stroke_color(curve.color))
    print(buf, "2.0 w\n")
    print(buf, curve.dash ? "[7 5] 0 d\n" : "[] 0 d\n")
    if !isempty(curve.x)
        print(buf, @sprintf("%.3f %.3f m\n", tx(curve.x[1]), ty(curve.y[1])))
        for j in 2:length(curve.x)
            print(buf, @sprintf("%.3f %.3f l\n", tx(curve.x[j]), ty(curve.y[j])))
        end
        print(buf, "S\n")
    end
    print(buf, "[] 0 d\n")
    return String(take!(buf))
end

function write_pdf(filename::String, content::String, W::Float64, H::Float64)
    objects = String[]
    push!(objects, "<< /Type /Catalog /Pages 2 0 R >>")
    push!(objects, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    push!(objects, @sprintf("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.0f %.0f] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>", W, H))
    push!(objects, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    push!(objects, "<< /Length $(ncodeunits(content)) >>\nstream\n$(content)\nendstream")

    open(filename, "w") do io
        write(io, "%PDF-1.4\n")
        offsets = Int[]
        for (i, obj) in enumerate(objects)
            push!(offsets, position(io))
            write(io, "$(i) 0 obj\n")
            write(io, obj)
            write(io, "\nendobj\n")
        end
        xref_pos = position(io)
        write(io, "xref\n")
        write(io, "0 $(length(objects)+1)\n")
        write(io, "0000000000 65535 f \n")
        for off in offsets
            write(io, @sprintf("%010d 00000 n \n", off))
        end
        write(io, "trailer\n")
        write(io, "<< /Size $(length(objects)+1) /Root 1 0 R >>\n")
        write(io, "startxref\n")
        write(io, "$(xref_pos)\n")
        write(io, "%%EOF\n")
    end
end

function write_pdf_plot(filename::String, curves::Vector{Curve})
    xmin, xmax = 0.0, 0.5
    ymin, ymax = 0.85, 3.15

    W, H = 900.0, 620.0
    left, right = 95.0, 35.0
    bottom, top = 85.0, 55.0
    plotW = W - left - right
    plotH = H - bottom - top

    function tx(x)
        return left + (x - xmin)/(xmax - xmin)*plotW
    end
    function ty(y)
        return bottom + (y - ymin)/(ymax - ymin)*plotH
    end

    buf = IOBuffer()
    print(buf, "1 1 1 rg 0 0 $(W) $(H) re f\n")

    if SHOW_TITLE
        print(buf, set_fill_color((0.0,0.0,0.0)))
        print(buf, pdf_text(W/2, H-30, 18, "Running interaction dimension on the square torus", align=:center))
    end

    # Grid and axes.
    print(buf, "0.85 0.85 0.85 RG 0.5 w\n")
    for x in 0.0:0.1:0.5
        print(buf, line_cmd(tx(x), bottom, tx(x), bottom+plotH))
    end
    for y in 1.0:0.5:3.0
        print(buf, line_cmd(left, ty(y), left+plotW, ty(y)))
    end

    print(buf, "0 0 0 RG 1.2 w\n")
    print(buf, line_cmd(left, bottom, left+plotW, bottom))
    print(buf, line_cmd(left, bottom, left, bottom+plotH))
    print(buf, line_cmd(left+plotW, bottom, left+plotW, bottom+plotH))
    print(buf, line_cmd(left, bottom+plotH, left+plotW, bottom+plotH))

    # Tick marks and labels.
    print(buf, set_fill_color((0.0,0.0,0.0)))
    print(buf, "0 0 0 RG 1 w\n")
    for x in 0.0:0.1:0.5
        px = tx(x)
        print(buf, line_cmd(px, bottom, px, bottom-5))
        print(buf, pdf_text(px, bottom-24, 11, @sprintf("%.1f", x), align=:center))
    end
    for y in 1.0:0.5:3.0
        py = ty(y)
        print(buf, line_cmd(left-5, py, left, py))
        label = (abs(y-round(y)) < 1e-10) ? @sprintf("%.0f", y) : @sprintf("%.1f", y)
        print(buf, pdf_text(left-12, py-4, 11, label, align=:right))
    end

    # Axis labels.
    print(buf, pdf_text(left + plotW/2, 28.0, 14, "r/L", align=:center))
    if SHOW_Y_AXIS_LABEL
        print(buf, pdf_text(12.0, bottom + plotH/2, 14, "D_delta", align=:left))
    end

    # Curves.
    for curve in curves
        print(buf, polyline_cmd(curve, xmin, xmax, ymin, ymax, left, bottom, plotW, plotH))
    end

    # Legend.
    legend_x = left + 20.0
    legend_y = bottom + plotH - 25.0
    dy = 18.0
    for (i, curve) in enumerate(curves)
        y = legend_y - (i-1)*dy
        print(buf, set_stroke_color(curve.color))
        print(buf, "2 w\n")
        print(buf, curve.dash ? "[7 5] 0 d\n" : "[] 0 d\n")
        print(buf, line_cmd(legend_x, y+4, legend_x+32, y+4))
        print(buf, "[] 0 d\n")
        print(buf, set_fill_color((0.0,0.0,0.0)))
        print(buf, pdf_text(legend_x+40, y, 10, curve.label, align=:left))
    end

    write_pdf(filename, String(take!(buf)), W, H)
end

# ------------------------------ output ---------------------------------

function write_csv(filename::String, xis::Vector{Float64}, Dtorus::Vector{Float64}, Dinf::Vector{Float64})
    open(filename, "w") do io
        println(io, "xi,D_torus,D_infinite")
        for j in eachindex(xis)
            @printf(io, "%.16e,%.16e,%.16e\n", xis[j], Dtorus[j], Dinf[j])
        end
    end
end

function main()
    mkpath(OUT_DIR)
    xis = collect(range(XI_MIN, XI_MAX, length=N_XI))

    colors = [
        (0.10, 0.20, 0.70),
        (0.75, 0.10, 0.10),
        (0.10, 0.55, 0.20)
    ]

    curves = Curve[]
    report_lines = String[]
    push!(report_lines, "Running interaction dimension convergence report")
    timestamp = Dates.format(Dates.now(), Dates.DateFormat("yyyy-mm-dd HH:MM:SS"))
    push!(report_lines, "Generated: $timestamp")
    push!(report_lines, "Formula: D_T2(xi;a) = 3 - xi*dU_T2/dxi/U_T2")
    push!(report_lines, "Reference: D_inf(s) = 3 - 2*s^2/((1+s^2)*log(1+s^2)), s=xi/a")
    push!(report_lines, @sprintf("potential-tail tolerance = %.3e", TAIL_TOL))
    push!(report_lines, @sprintf("dimension convergence tolerance = %.3e", DIM_CONV_TOL))
    push!(report_lines, "")

    for (ia, a) in enumerate(A_VALUES)
        println("Processing a = ", a)
        N, tail, derr = choose_cutoff_for_dimension(xis, a)
        Dtorus = torus_dimension_axial(xis, a, N)
        Dinf = [infinite_dimension_reference_xi(xi, a) for xi in xis]

        csv_file = joinpath(OUT_DIR, "running_dimension_a_$(sanitize_a(a)).csv")
        write_csv(csv_file, xis, Dtorus, Dinf)

        col = colors[ia]
        push!(curves, Curve(xis, Dtorus, @sprintf("torus a=%.2f", a), col, false))
        push!(curves, Curve(xis, Dinf, @sprintf("plane a=%.2f", a), col, true))

        @printf("  cutoff N = %d, tail = %.3e, max |D_2N - D_N| = %.3e\n", N, tail, derr)
        push!(report_lines, @sprintf("a = %.6f", a))
        push!(report_lines, @sprintf("  cutoff N                  = %d", N))
        push!(report_lines, @sprintf("  tail estimate             = %.8e", tail))
        push!(report_lines, @sprintf("  max |D_2N - D_N|          = %.8e", derr))
        push!(report_lines, "  csv                       = $(basename(csv_file))")
        push!(report_lines, "")
    end

    pdf_file = joinpath(OUT_DIR, "running_dimension_torus.pdf")
    write_pdf_plot(pdf_file, curves)

    report_file = joinpath(OUT_DIR, "running_dimension_convergence_report.txt")
    open(report_file, "w") do io
        for line in report_lines
            println(io, line)
        end
        println(io, "Files written:")
        println(io, "  running_dimension_torus.pdf")
        for a in A_VALUES
            println(io, "  running_dimension_a_$(sanitize_a(a)).csv")
        end
    end

    println("Done. Output directory: ", OUT_DIR)
    println("PDF figure: ", pdf_file)
    println("Convergence report: ", report_file)
end

main()
