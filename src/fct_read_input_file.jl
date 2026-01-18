function read_input_file(filepath)
    f = open(filepath, "r")

    # Skip the first three lines
    for i = 1:3
        readline(f)
    end

    # Read the remaining lines
    lines = []
    while !eof(f)
        line = readline(f)
        if occursin("#", line) && !isempty(line)
            line = strip(split(line, "#")[1])
        end
        if !startswith(line, "#") && !isempty(line)
            push!(lines, split(line, ":")[2])
        end
    end

    # Close the file
    close(f)

    # Return lines
    return(strip(lowercase(lines[1])), 
           parse(Int, lines[2]), 
           parse(Int, lines[3]), 
           parse(Int, lines[4]), 
           parse(Int, lines[5]),
           parse(Int, lines[6]),
           parse(Int, lines[7]),
           parse(Int, lines[8]),
           parse(Int, lines[9]),
           parse.(Int, split(lines[10],",")),
           strip(lines[11]),
           strip(lines[12]),
           strip(lines[13]),
           strip(lines[14]),
           strip(lines[15]),
           strip(lines[16]),
           strip(lines[17]),
           strip(lines[18]),
           strip(lines[19]),
           parse(Int, lines[20]))
end


