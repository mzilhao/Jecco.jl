
using Dates

function startup()

    date        = string(now())
    host        = gethostname()
    user        = ENV["USER"]
    num_threads = ENV["JULIA_NUM_THREADS"]

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
