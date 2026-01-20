using FileIO,JLD2
using CairoMakie
using CairoMakie: translate!
using Colors

theme = Theme(
    font = "Poppins Regular" ,
    Axis = (;
        xticklabelsize = 16, 
        xticklabelfont = "Poppins Regular",
        yticklabelsize = 16, 
        yticklabelfont = "Poppins Regular",
        xlabelsize = 20, 
        xlabelfont = "Poppins Regular",
        ylabelsize = 20, 
        ylabelfont = "Poppins Regular",
        titlefont = "Poppins Medium", 
        titlesize = 20
        ), 
    Label = (; font = "Poppins Regular"),
    Colorbar = (;labelfont = "Poppins Regular", labelsize = 20)
)
set_theme!(theme)

path = "FigureSources/figure_S2.jld2"

vars = load(path)
for key in ["baseM","maxL","freq","w","μ_range","θ","δ_list","delta_range","δ_num","h","numPoints_1","numPoints_2","numPoints_L","neighbours","switch_guess_points","solTolerance","uniquenessTolerance","solutions","numsols",]
    @eval $(Symbol(key)) = $(vars[key])
end

f = Figure(;size = (1000,1000))
ga = f[1,1] = GridLayout()
gb = f[2,1] = GridLayout()

t = 3
before = Axis(ga[1,1])
axA = Axis(ga[1,2:t]; xlabel = "Shape μ", ylabel = "Phenotypic Spread δ", aspect = 1)
limits!(axA, [0,1],[0.5,15.5])
axA.xticks = 0:0.1:1
axA.yticks = 1:1:15
cf = contourf!(axA,μ_range, delta_range,numsols;levels = [0;1.5;3;5;7;9;11;13;15;17;19;21;23;25;27;29;31],colormap = :turbo)
cbar = Colorbar(ga[1,t+1], cf;label = "No. of Solutions", ticks = sort(unique(numsols)))#, height = @lift Fixed($(pixelarea(axA.scene)).widths[1]))
after = Axis(ga[1,t+2])
hidedecorations!.((before,after))
hidespines!.((before,after))

gridaxes = [
    Axis(gb[i,j]; xlabel= (i==4 && j==2) ? "Day of Lunar Cycle" : (""), ylabel = (i==2 && j==1) ? "Fraction Emerging" : "") 
    for i in 1:4, j in 1:4
    ]

for i in 1:4
    for j in 1:4
        gridaxes[i,j].xticks = [-5;0;5;10;15;20;25]
        if j ∈ [1,2]
            ylims!(gridaxes[i,j], [nothing,0.79])
        elseif i ∈ [1,2]
            ylims!(gridaxes[i,j], [nothing,0.08])
        elseif i ∈ [3,4]
            ylims!(gridaxes[i,j], [nothing,0.22])
            gridaxes[i,j].yticks = [0;0.1;0.2]
        end
    end
end

hidexdecorations!.(gridaxes[1:end-1,:]; ticks = false, grid = false, minorgrid = false, minorticks = false)
# hideydecorations!.(gridaxes[:,2]; ticks = false, grid = false, minorgrid = false, minorticks = false)
# hideydecorations!.(gridaxes[:,4]; ticks = false, grid = false, minorgrid = false, minorticks = false)

μs = [30; 50; 85; 99]
δs = [2; 5; 9; 14]
for (xpos,μ_idx) in enumerate(μs)
    for (ypos, δ_idx) in enumerate(reverse(δs))
        barplot!(gridaxes[ypos,xpos], -5:24, solutions[:,μ_idx,δ_idx])
    end
end

for x in μs
    for y in δs
        scatter!(axA,x/100,y; color = :white, markersize = 20)
        scatter!(axA,x/100,y; color = :black, markersize = 10)
    end
end

Label(ga[1, 2, TopLeft()], "A",fontsize = 30,font = :bold,padding = (-20, 0, 0, 30), halign = :left, tellheight = false)
Label(gb[1, 1, TopLeft()], "B",fontsize = 30,font = :bold,padding = (0, 0, 0, 0), halign = :left, tellheight = false)

rowsize!(f.layout, 1, 300)
colgap!(gb,10)
rowgap!(gb,15)


translate!(gridaxes[2,1].yaxis.elements[:labeltext], 0, -75)
translate!(gridaxes[4,2].xaxis.elements[:labeltext], 125, 0)

try
    Base.Filesystem.mkdir("Figures")
catch
end
save("Figures/figure_S2.pdf", f)
f