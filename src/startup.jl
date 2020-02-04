
using Dates

function startup()

    date        = string(now())
    host        = gethostname()
    user        = ENV["USER"]

    num_threads = try
        ENV["JULIA_NUM_THREADS"]
    catch e
        if isa(e, KeyError)
            # no JULIA_NUM_THREADS defined; set num_threads to 1
            1
        else
            # unknown error
            error(e)
        end
    end

    println("-------------------------------------------------------------")
    println("")

    println(raw"   ^..^        // ")
    println(raw"   /_/\\========        Jecco")
    println(raw"      //\\   //\\  ")
    println(raw"     //  \\ //  \\ ")
    println("")

    println("-------------------------------------------------------------")
    println("")

    println("Run date:          $date")
    println("Run host:          $host")
    println("username:          $user")
    println("")
    println("This process contains $num_threads threads")
    println("")

end
