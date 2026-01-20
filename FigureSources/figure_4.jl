using JLD2
using CairoMakie
using Colors

theme = Theme(
    font = "Poppins Semibold" ,
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

path = "FigureSources/figure_4.jld2"

vars = load(path)
for key in ["numsols", "μ_range", "δ_list", "delta_list", "baseM","maxL","freq","w","θ","h","numPoints_1","numPoints_L","neighbours","solTolerance","uniquenessTolerance"]
    @eval $(Symbol(key)) = $(vars[key])
end

f = Figure(size = (1500, 700))
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()

PayoffFunction(x,h) = x*(1+h)/(1+h*x)
ax1 = Axis(ga[1,1], xlabel = "Type frequency x", ylabel = "Reproductive Benefit")
x = range(0,1,100)

lines!(ax1, x, PayoffFunction.(x,0), linewidth = 10, color = :black)
for h in 1:2:9
    lines!(ax1, x, PayoffFunction.(x,h),linewidth = 9, color = RGBf(0.8,0.8,0.8))
end
lines!(ax1, x, PayoffFunction.(x,10), linewidth = 10, color = :black)    

text!(0.56,0.34; text = "h=0\n(Fig. 3)", fontsize = 26, align=(:center,:bottom), color = (:black, 0.5))
text!(0.14,0.81; text = "h=10", fontsize = 26)

fig = Figure()
ax2 = Axis(gb[1,1],xlabel = "Shape μ", ylabel = "Phenotypic Spread δ", title = "h=10")
limits!(ax2, [0,1], [0.5,15.5])
ax2.xticks = 0:0.1:1
ax2.yticks = 1:1:15
hm = contourf!(ax2,μ_range,delta_list,numsols; levels = [0;1.5;3;5;7;9;11;13;15;17;19;21;23;25;27;29;31],colormap=:turbo)
Colorbar(gb[1,2],hm;label = "No. of Solutions", ticks = sort(unique(numsols)))

Label(ga[1, 1, TopLeft()], "A",fontsize = 30,font = :bold,padding = (0, 50, 0, 0),        halign = :center)
Label(gb[1, 1, TopLeft()], "B",fontsize = 30,font = :bold,padding = (0, 45, 0, 0),        halign = :center)

colsize!(f.layout, 1, Auto(0.9))

## Drawing the border of the regions
numsols_border = deepcopy(numsols)
numsols_border_limited = @view numsols_border[:, 7:end]
numsols_border_limited[numsols_border_limited .< 3] .= 4

border_kwargs = (; color = :white, linestyle = (:dash, 2.5), linewidth = 7, linecap = :round, joinstyle = :round)
cf2 = contour!(ax2, μ_range, delta_list, numsols_border; border_kwargs..., levels = [3])
cf3 = contour!(ax2, μ_range, delta_list, numsols; border_kwargs..., levels = [1.5])
μ_min_1 = μ_range[findfirst(==(1), numsols[:, end])]
μ_min_2 = μ_range[findfirst(==(2), numsols[:, 1])]

lines!(ax2, [(μ_min_1, 15.45), (0.998, 15.45), (0.998, 10.5)]; border_kwargs...)
lines!(ax2, [(μ_min_2, 0.55), (0.998, 0.55), (0.998, 6.5)]; border_kwargs...)

# Labelling the regions
text!(ax2, 0.65, 3; text = "Lunar", color = :white, fontsize = 26, font = "Poppins Semibold")
text!(ax2, 0.74, 13.8; text = "Arrhythmic", color = :white, fontsize = 26, font = "Poppins Semibold")

try
    Base.Filesystem.mkdir("Figures")
catch
end

save("Figures/figure_4.pdf",f)
f