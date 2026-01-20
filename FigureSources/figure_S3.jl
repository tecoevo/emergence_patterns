using JLD2
using CairoMakie
using Colors
using LazySets: Interval

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
    Colorbar = (;labelfont = "Poppins Regular", labelsize = 26))
    set_theme!(theme)

path = "FigureSources/figure_S3.jld2"
vars = load(path)
for key in ["baseM","maxL","freq","μ_begin","μ_end","μ_num","θ","δ_list","δ_num","h","numPoints_1","numPoints_2","numPoints_L","neighbours","solTolerance","uniquenessTolerance","μ_range","δ_range","delta_list", "numsols_0_1", "numsols_0_4", "numsols_0_6"]
    @eval $(Symbol(key)) = $(vars[key])
end


cp_levels = [0;1.5;3;5;7;9;11;13;15;17;19;21;23;25;27;29;31]
cbar_levels = [1;2;4;6;8;10;12;14;16;18;20;22;24;26;28;30]

f = Figure(size = (1500,500))
ax1 = Axis(f[1,1],xlabel = "Shape μ", ylabel = "Phenotypic Spread δ", title = "𝗔. w = 0.1")
limits!(ax1, [0,1],[0.5,15.5])
ax1.xticks = 0:0.1:1
ax1.yticks = 1:1:15
# ct = contourf!(ax,μ_range, delta_list,numsols_0_1;levels = cp_levels,colormap = :turbo)
hm1 = heatmap!(ax1,μ_range, δ_range , numsols_0_1; colormap = :turbo)

ax2 = Axis(f[1,2],xlabel = "Shape μ", title = "𝗕. w = 0.4")
limits!(ax2, [0,1],[0.5,15.5])
ax2.xticks = 0:0.1:1
ax2.yticks = 1:1:15
hm2 = heatmap!(ax2,μ_range, δ_range , numsols_0_4; colormap = :turbo)    

ax3 = Axis(f[1,3],xlabel = "Shape μ", title = "𝗖. w = 0.6")
limits!(ax3, [0,1],[0.5,15.5])
ax3.xticks = 0:0.1:1
ax3.yticks = 1:1:15
hm3 = heatmap!(ax3,μ_range, δ_range , numsols_0_6; colormap = :turbo)    

Colorbar(f[1,4]; colormap = cgrad(:turbo, cbar_levels./30, categorical = true), limits = (0.5,31), label = "No. of Solutions", ticks = cbar_levels)

try
    Base.Filesystem.mkdir("Figures")
catch
end
save("Figures/figure_S3.pdf", f)
f