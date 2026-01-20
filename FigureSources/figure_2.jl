using CairoMakie
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
        titlesize = 22
        ), 
    Label = (; font = "Poppins Regular"),
    Colorbar = (;labelfont = "Poppins Regular", labelsize = 18))
set_theme!(theme)

try
    Base.Filesystem.mkdir("Figures")
catch
end

begin #barplot of just the payoffs
    function cosinePayoffVector(w,θ,M,freq=1) 
        A = 1 .+ w * cos.(2*freq*pi/M*(collect(1:M) .- θ))
        A[A .< 0.] .= 0.
        return A
    end
    payoff = cosinePayoffVector(1, 8.25, 30, 2)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel = "Day of Lunar Cycle", ylabel = "Payoff")
    ylims!(ax,[0,nothing])
    hidexdecorations!(ax,label=false,ticklabels=false,ticks=false)
    hideydecorations!(ax,label=false,ticklabels=false,ticks=false)
    barplot!(ax, payoff,color=Colors.HSV(44,77,100))
    save("Figures/payoffs.svg",fig)
    display(fig)
end

begin
    fig2 = Figure()
    ax2 = Axis(fig2[1,1], xlabel = "Type frequency x", ylabel = "Reproductive Benefit")
    x = range(0,1,100)
    lines!(ax2, x, x, linewidth = 7.5)
    save("Figures/reproductive_benefit.svg",fig2)
    display(fig2)
end

begin
    fig3 = Figure(;size = (1000,450))
    ax3 = Axis(fig3[1,1]; xlabel = "Type", ylabel = "Proportion of Offspring")
    ax3.xticks = [1;3;5;7;9]
    ax4 = Axis(fig3[1,2]; xlabel = "Type", ylabel = "Proportion of Offspring")
    ax4.xticks = [1;3;5;7;9]
    y1 = [0;0;1;3;9;4;1;0;0]
    y2 = [0;0;6;7.7;9;8;5;0;0]
    hideydecorations!.((ax3,ax4); label = false)
    hidexdecorations!.((ax3,ax4); label = false, ticks = false, ticklabels= false)
    ylims!(ax3,[0,nothing])
    ylims!(ax4,[0,nothing])
    hidespines!.((ax3,ax4), :r)
    hidespines!.((ax3,ax4), :t)
    barplot!(ax3,y1)
    barplot!(ax4,y2)
    save("Figures/phenotypic_distribution.svg",fig3)
    display(fig3)
end