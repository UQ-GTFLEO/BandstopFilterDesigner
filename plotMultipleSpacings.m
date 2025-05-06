spacing = 1:0.5:5;
for space = spacing
    generate_required_impedance_curve(0:0.1:40, space);
end