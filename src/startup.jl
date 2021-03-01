
using Dates

function startup(outdir::String; remove_existing::Bool=false)
    # if no outdir specified, error out
    if outdir == ""
        error("Output folder cannot be empty string.")
    end

    # create output folder if it doesn't exist already
    if isdir(outdir) && remove_existing
        rm(outdir, recursive=true)
    end
    if !isdir(outdir)
        mkdir(outdir)
    end

    # copy execution script to outdir folder
    scriptname = basename(Base.source_path())
    dst = joinpath(outdir, scriptname)
    try
        cp(Base.source_path(), dst)
    catch e
        if isa(e, ArgumentError) # file already exist
            nothing
        else
            throw(e)
        end
    end

    timenow = now()
    datef   = Dates.format(timenow, "e, dd u yyyy (HH:MM:SS)")

    date    = string(datef)
    host    = gethostname()
    user    = try
        ENV["USER"]
    catch e
        if isa(e, KeyError) # probably on windows
            splitdir(homedir())[end]
        else
            throw(e)
        end
    end

    num_threads = try
        ENV["JULIA_NUM_THREADS"]
    catch e
        if isa(e, KeyError)
            # no JULIA_NUM_THREADS defined; set num_threads to 1
            1
        else
            # unknown error
            throw(e)
        end
    end

    println("-------------------------------------------------------------")

    println(raw"                                     _      ")
    println(raw"                                 .:'/   _..._  ")
    println(raw"           Jecco                // ( ```-.._.'  ")
    println(raw"                                \| /    0\___    ")
    println(raw"                                |     0      \    ")
    println(raw"       Julia Einstein           |            /  ")
    println(raw"                                \_       .--'  ")
    println(raw"    Characteristic Code         (_'---'`)      ")
    println(raw"                                / `'---`|  ")
    println(raw"                              ,'        |  ")
    println(raw"              .            .'`          |  ")
    println(raw"              )\       _.-'             ;  ")
    println(raw"             / |    .'`   _            /  ")
    println(raw"           /` /   .'       '.        , |  ")
    println(raw"          /  /   /           \   ;   | |  ")
    println(raw"          |  \  |            |  .|   | |  ")
    println(raw"           \  ``|           /.-' |   | |  ")
    println(raw"            '-..-\       _.;.._  |   |.;-.  ")
    println(raw"                  \    <`.._  )) |  .;-. ))  ")
    println(raw"                  (__.  `  ))-'  \_    ))'   ")
    println(raw"                      `'--``       `````  ")

    println("-------------------------------------------------------------")
    println("")

    println("Run date:          $date")
    println("Run host:          $host")
    println("username:          $user")
    println("")
    println("This process contains $num_threads threads")
    println("")

end
