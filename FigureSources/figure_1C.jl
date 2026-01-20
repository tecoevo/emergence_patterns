using CSV
using CircularArrays
using CairoMakie

theme = Theme(
    font = "Poppins Regular" ,
    Axis = (;
        xticklabelsize = 20, 
        xticklabelfont = "Poppins Regular",
        yticklabelsize = 20, 
        yticklabelfont = "Poppins Regular",
        xlabelsize = 26, 
        xlabelfont = "Poppins Regular",
        ylabelsize = 26, 
        ylabelfont = "Poppins Regular",
        titlefont = "Poppins Medium", 
        titlesize = 36
        ), 
    Label = (;font = "Poppins Regular"),
    Colorbar = (;labelsize = 28))
set_theme!(theme)

file = CSV.File("FigureSources/figure_1C.csv")
days = CircularArray(file.days)
lunar = file.lunar
semilunar = file.semilunar
arrhythmic = file.arrhythmic

g = Figure(; size = (450, 700))
f = g[1,1] = GridLayout()
kwargs = (; xminorgridvisible = true, xminorgridcolor = RGBAf(0, 0, 0, 0.1), xticks = 0:10:20, xminorticks = [-5, 5, 15, 25], xminorticksvisible = true, xtickwidth = 2)
ax1 = Axis(f[1,1]; kwargs...)
ax2 = Axis(f[2,1]; kwargs...)
ax3 = Axis(f[3,1]; kwargs...)
ax4 = Axis(f[4,1]; xlabel = "Day of Lunar Cycle", xminorgridvisible = true, kwargs...)
hidexdecorations!.((ax1, ax2, ax3),ticks=false, grid = false, minorgrid = false, minorticks = false)
ylab = Label(f[1:4, 0], "Proportion Emerging", rotation = pi/2, fontsize = 26, padding = (-6, 10, 0.0f0, 0.0f0))

barplot!(ax1, days, lunar, color = :steelblue3)
barplot!(ax2, days[16:45], lunar, color = :magenta4)
barplot!(ax3, days[-2:27], semilunar, color = :darkgoldenrod1)
barplot!(ax4, days, arrhythmic, color = :lightgreen)

ylims!(ax3, nothing, 0.23)
ylims!(ax4, nothing, 0.055)
hideydecorations!.((ax1, ax2, ax3, ax4); label = false)

Box(f[1,2], color = :gray90, width = 45)
Box(f[2,2], color = :gray90, width = 45)
Box(f[3,2], color = :gray90, width = 45)
Box(f[4,2], color = :gray90, width = 45)
Label(f[1,2], "Lunar\nFull Moon", rotation = pi/2, tellheight = false)
Label(f[2,2], "Lunar\nNew Moon", rotation = pi/2, tellheight = false)
Label(f[3,2], "Semilunar", rotation = pi/2, tellheight = false)
Label(f[4,2], "Arrhythmic", rotation = pi/2, tellheight = false)

colgap!(f,2,0)
colgap!(f,1,0)
rowgap!(f,1,15)
rowgap!(f,2,15)
rowgap!(f,3,15)

Label(
    g[1, 1, TopLeft()], "C",
    fontsize = 28,
    font = "Poppins Bold",
    padding = (4, 0, -20, 0),
    halign = :right,
    valign = :bottom,
    width = 0 
)

display(g)

try
    Base.Filesystem.mkdir("Figures")
catch
end
save("Figures/figure_1C.pdf", g) # needs to be combined with the svg for A, B