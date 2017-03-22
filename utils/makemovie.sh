#!/bin/bash
function usage
{
    echo "usage: makemovie.sh [-i files] [-o file] [-f fps] [-h]"
}

INPUT="%06d.jpg"
OUTPUT=output.mp4
FPS=30
while [ "$1" != "" ]; do
    case $1 in
        -i | --input )          shift
                                INPUT=$1
                                ;;
        -o | --output )         shift
                                OUTPUT=$1
                                ;;
        -f | --fps )            shift
                                FPS=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

ffmpeg -framerate $FPS -i $INPUT -c:v libx264 -vf format=yuv420p $OUTPUT 
