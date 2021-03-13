using Meshes
using Printf


function check_file(filename)
    overwrite = true

    if isfile(filename)
        fname = basename(filename)
        print("File " * fname * " already exists in the chosen directory. ")
        print("Do you wish to overwrite it? [Y/n] ")
        choice = readline()

        if lowercase(first(choice)) == "n"
            overwrite = false
        end
    end

    overwrite
end

const MESH_FORMATS = ["Compact", "GMSH"]

function save_mesh(g::CartesianGrid, filename::AbstractString; format="Compact")

    if !(format in MESH_FORMATS)
        print("Mesh was not saved: Format not supported")
    elseif !check_file(filename)
        print("Mesh was not saved: File already exists.")

    elseif format == "Compact"

        open(filename, "w") do f
            write(f, "CartesianGrid\n")

            minim = string.( minimum(g).coords )
            write(f, "Minimum: " * join(minim, " ") * '\n')

            maxim = string.( maximum(g).coords )
            write(f, "Maximum: " * join(minim, " ") * '\n')

            step = string.( g.spacing )
            write(f, "Spacing: " * join(step, " "))
        end

    else
        save_gmsh(g, filename, file_checked=true)
    end

    nothing
end

function save_mesh(m::Unstructuredmesh, filename::AbstractString; format="GMSH")
    if format != "GMSH"
        print("Mesh was not saved: Format not supported")
    elseif !check_file(filename)
        print("Mesh was not saved: File already exists.")
    else
        save_gmsh(m, filename, file_checked=true)
    end
    nothing
end

function save_gmsh(m::Mesh, filename::AbstractString; file_checked=false)
    nothing
end
