function parse_intraday_hours(s::AbstractString, max_hour::Int=48)
    s = strip(uppercase(s))

    if s == "N"
        return Int[]
    end

    if s == "ALL"
        return collect(1:max_hour)
    end

    hours = Int[]
    for token in split(s, ",")
        token = strip(token)
        isempty(token) && continue
        if occursin("-", token)
            parts = split(token, "-")
            length(parts) == 2 || error("Invalid range '$token'. Use 'a-b' format.")
            a, b = parse(Int, strip(parts[1])), parse(Int, strip(parts[2]))
            a <= b || error("Invalid range '$token': start must be ≤ end.")
            append!(hours, a:b)
        else
            push!(hours, parse(Int, token))
        end
    end
    return hours
end


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
    scenario_length = parse(Int, lines[2])
    return (
        # Data type - load, solar, wind
        strip(lowercase(lines[1])),
        # Scenario length
        scenario_length,
        # Scenario paths/number of scenarios
        parse(Int, lines[3]),
        #
        parse(Int, lines[4]),
        #
        parse(Int, lines[5]),
        #
        parse(Int, lines[6]),
        #
        parse(Int, lines[7]),
        #
        parse(Int, lines[8]),
        #
        parse(Int, lines[9]),
        #
        parse_intraday_hours(lines[10], scenario_length),
        # 
        strip(lines[11]),
        #
        strip(lines[12]),
        #
        strip(lines[13]),
        #
        strip(lines[14]),
        #
        strip(lines[15]),
        #
        strip(lines[16]),
        #
        strip(lines[17]),
        #
        strip(lines[18]),
        #
        strip(lines[19]),
        #
        parse(Int, lines[20]))
end


