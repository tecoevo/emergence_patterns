using JLD2
using CairoMakie
using Colors

theme = Theme(
    font = "Poppins Regular",
    Axis = (;
        xticklabelsize = 20, 
        xticklabelfont = "Poppins Regular",
        yticklabelsize = 20, 
        yticklabelfont = "Poppins Regular",
        xlabelsize = 28, 
        xlabelfont = "Poppins Regular",
        ylabelsize = 28, 
        ylabelfont = "Poppins Regular",
        titlefont = "Poppins Medium", 
        titlesize = 22
        ), 
    Label = (; font = "Poppins Regular"),
    Colorbar = (;labelfont = "Poppins Regular", labelsize = 26)
)
set_theme!(theme)

path = "FigureSources/figure_3.jld2"

vars = load(path)
for key in ["numsols","arrhythmic_solution","fullmoon_solution","newmoon_solution", "μ_range", "δ_list", "delta_list", "baseM","maxL","freq","w","θ","h","numPoints_1","numPoints_L","neighbours","solTolerance","uniquenessTolerance"]
    @eval $(Symbol(key)) = $(vars[key])
end

f = Figure(size = (1500, 700))
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()

ax1 = Axis(ga[1,1],xlabel = "Shape μ", ylabel = "Phenotypic Spread δ")
limits!(ax1, [0,1],[0.5,15.5])
ax1.xticks = 0:0.1:1
ax1.yticks = 1:1:15
cf = contourf!(ax1, μ_range, delta_list, numsols; levels = [0;1.5;3;5;7;9;11;13;15;17;19;21;23;25;27;29;31], colormap = :turbo)
# Colorbar(ga[1,2]; colormap = cgrad(:turbo, sort(unique(numsols))./30, categorical = true), limits = (0.5,31), label = "No. of Solutions", ticks = sort(unique(numsols)))
Colorbar(ga[1,2], cf; label = "No. of Solutions", ticks = sort(unique(numsols)))

ax2 = Axis(gb[1,1])
ax3 = Axis(gb[2,1], ylabel = "Fraction Emerging")
ax4 = Axis(gb[3,1], xlabel = "Day of Lunar Cycle")
ax3.yticklabelspace = 50. 
ax2.xticks = [-5;0;5;10;15;20;25]
ax3.xticks = [-5;0;5;10;15;20;25]
ax4.xticks = [-5;0;5;10;15;20;25]

hidexdecorations!(ax2; ticks = false, grid = false, minorgrid = false, minorticks = false)
hidexdecorations!(ax3; ticks = false, grid = false, minorgrid = false, minorticks = false)

barplot!(ax2, -5:24, circshift(arrhythmic_solution,-2), color = :lightgreen)
barplot!(ax3, -5:24, circshift(fullmoon_solution,-2), color = RGB(0.32421875, 0.54296875, 0.73046875))
barplot!(ax4, -5:24, circshift(newmoon_solution,-2), color = RGB(0.453,0.078, 0.48))

Box(gb[1,2], color = :gray90)
Box(gb[2,2], color = :gray90, width = 45)
Box(gb[3,2], color = :gray90)
Label(gb[1,2], "Arrhythmic", rotation = pi/2, tellheight = false)
Label(gb[2,2], "Lunar\nFull Moon", rotation = pi/2, tellheight = false)
Label(gb[3,2], "Lunar\nNew Moon", rotation = pi/2, tellheight = false)
colgap!(gb,1,0)

xspace = maximum(tight_xticklabel_spacing!,[ax1,ax4])
ax1.xticklabelspace = xspace
ax4.xticklabelspace = xspace

Label(ga[1, 1, TopLeft()], "A", fontsize = 30,font = :bold,padding = (0, 50, 5, 0), halign = :right)
Label(gb[1, 1, TopLeft()], "B", fontsize = 30,font = :bold,padding = (0, 80, 5, 0), halign = :right)

colsize!(f.layout, 2, Auto(0.7))

## Drawing the border of the regions
numsols_border = deepcopy(numsols)
numsols_border_limited = @view numsols_border[:, 7:end]
numsols_border_limited[numsols_border_limited .< 3] .= 4

border_kwargs = (; color = :white, linestyle = (:dash, 2.5), linewidth = 7, linecap = :round, joinstyle = :round)
cf2 = contour!(ax1, μ_range, delta_list, numsols_border; border_kwargs..., levels = [3])
cf3 = contour!(ax1, μ_range, delta_list, numsols; border_kwargs..., levels = [1.5])
μ_min_1 = μ_range[findfirst(==(1), numsols[:, end])]
μ_min_2 = μ_range[findfirst(==(2), numsols[:, 1])]

lines!(ax1, [(μ_min_1, 15.45), (0.998, 15.45), (0.998, 10.5)]; border_kwargs...)
lines!(ax1, [(μ_min_2, 0.55), (0.998, 0.55), (0.998, 6.5)]; border_kwargs...)

# Labelling the regions
text!(ax1, 0.79, 3; text = "Lunar", color = :white, fontsize = 26, font = "Poppins Semibold")
text!(ax1, 0.815, 14.5; text = "Arrhythmic", color = :white, fontsize = 26, rotation = -π/4-0.1, font = "Poppins Semibold")

try
    Base.Filesystem.mkdir("Figures")
catch
end

save("Figures/figure_3.pdf", f)
f