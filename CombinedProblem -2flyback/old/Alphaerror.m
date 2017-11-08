function error = Alphaerror(lift_search,v,V,scattered,Alpha_desired)

Alpha = scattered.AoA(v,V,lift_search);

error = abs(Alpha_desired - Alpha);

end
