using JLD2
using CairoMakie
using CairoMakie: translate!
using Colors

theme = Theme(
    font = "Poppins Regular" ,
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

path = "FigureSources/figure_5.jld2"

vars = load(path)
for key in ["numsols","arrhythmic_solution","semilunar_solution","fullmoon_solution","newmoon_solution", "μ_range", "δ_list", "δ_range", "delta_range", "baseM","maxL","freq","w","θ","h","numPoints_1","numPoints_L","neighbours","solTolerance","uniquenessTolerance"]
    @eval $(Symbol(key)) = $(vars[key])
end

f = Figure(size = (1500, 700))
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()

ax1 = Axis(ga[1,1],xlabel = "Shape μ", ylabel = "Phenotypic Spread δ")
limits!(ax1, [0,1],[0.5,15.5])
ax1.xticks = 0:0.1:1
ax1.yticks = 1:1:15
levels = [0.4;1.7;2.5;3.7;8;12;15;18;21;26;29;32;35;38;41.2;42.3;43.2;44.5;46]
ticks = (1.75:2.5:44.25,string.([1;2;3;6;10;14;16;20;24;28;30;34;36;40;42;43;44;45]))
hm = contourf!(ax1,μ_range,delta_range,numsols; levels = levels,colormap=:turbo)
Colorbar(ga[1,2]; colormap = cgrad(:turbo, 18, categorical = true), limits = (0.5,45.5), label = "No. of Solutions", ticks = ticks)

ax2 = Axis(gb[1,1])
ax3 = Axis(gb[2,1], ylabel = "Fraction Emerging")
ax4 = Axis(gb[3,1])
ax5 = Axis(gb[4,1], xlabel = "Day of Lunar Cycle")
ax3.yticklabelspace = 50. 
ax2.yticks = [0;0.02]
ax3.yticks = [0;0.1;0.2]
ax5.xticks = [-5;0;10;20;25]

hidexdecorations!(ax2; ticks = false, grid = false, minorgrid = false, minorticks = false)
hidexdecorations!(ax3; ticks = false, grid = false, minorgrid = false, minorticks = false)
hidexdecorations!(ax4; ticks = false, grid = false, minorgrid = false, minorticks = false)

barplot!(ax2, -5:24, circshift(arrhythmic_solution,-2), color = RGB(0.64453125, 0.921875, 0.59765625))
barplot!(ax3, -5:24, circshift(semilunar_solution,-2), color = RGB(0.9453125, 0.65234375, 0.23046875))
barplot!(ax4, -5:24, circshift(fullmoon_solution,-2), color = RGB(0.33203125, 0.52734375, 0.72265625))
barplot!(ax5, -5:24, circshift(newmoon_solution,-2), color = RGB(0.453125, 0.078125, 0.48046875))

Box(gb[1,2], color = :gray90)
Box(gb[2,2], color = :gray90)
Box(gb[3,2], color = :gray90, width = 45)
Box(gb[4,2], color = :gray90)
Label(gb[1,2], "Arrhythmic", rotation = pi/2, tellheight = false)
Label(gb[2,2], "Semilunar", rotation = pi/2, tellheight = false)
Label(gb[3,2], "Lunar\nFull Moon", rotation = pi/2, tellheight = false)
Label(gb[4,2], "Lunar\nNew Moon", rotation = pi/2, tellheight = false)
colgap!(gb,1,0)
rowgap!(gb,15)

xspace = maximum(tight_xticklabel_spacing!,[ax1,ax4])
ax1.xticklabelspace = xspace
ax4.xticklabelspace = xspace


Label(ga[1, 1, TopLeft()], "A",fontsize = 30,font = :bold,padding = (0, 50, 5, 0),        halign = :right)
Label(gb[1, 1, TopLeft()], "B",fontsize = 30,font = :bold,padding = (0, 80, 5, 0),        halign = :right)

colsize!(f.layout, 2, Auto(0.7))
translate!(ax3.yaxis.elements[:labeltext], 0, -75) # Anything that has field Transformation can be transformed by translate!, rotate! or scale!

## Drawing the border of the regions
numsols_border = deepcopy(numsols)
numsols_border_limited = @view numsols_border[:, 7:end]
numsols_border_limited[numsols_border_limited .< 4] .= 4

border_kwargs = (; color = :white, linestyle = (:dash, 2.5), linewidth = 7, linecap = :round, joinstyle = :round)
cf2 = contour!(ax1, μ_range, delta_range, numsols_border; border_kwargs..., levels = [3.7])
cf3 = contour!(ax1, μ_range, delta_range, numsols; border_kwargs..., levels = [2.5])
μ_min_1 = μ_range[findfirst(==(1), numsols[:, end])]
μ_min_2 = μ_range[findfirst(==(3), numsols[:, 1])]

lines!(ax1, [(μ_min_1, 15.45), (0.998, 15.45), (0.998, 10.5)]; border_kwargs...)
lines!(ax1, [(μ_min_2, 0.55), (0.998, 0.55), (0.998, 6.)]; border_kwargs...)

# Labelling the regions
text!(ax1, 0.84, 3.65; text = "Lunar\nor\nSemilunar", color = :white, fontsize = 26, font = "Poppins Semibold", align = (:center, :center))
text!(ax1, 0.815, 14.6; text = "Arrhythmic", color = :white, fontsize = 26, font = "Poppins Semibold", rotation = -π/4-0.1)

try
    Base.Filesystem.mkdir("Figures")
catch
end
save("Figures/figure_5.pdf", f)
f