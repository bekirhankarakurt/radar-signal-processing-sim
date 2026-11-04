# Terminal demo rendered with charmbracelet/vhs:  vhs docs/demo.tape
# Expects a Release build in ./build (cmake -S . -B build && cmake --build build -j).
Output docs/demo.gif

Set Shell "bash"
Set FontSize 15
Set Width 1100
Set Height 640
Set Padding 16
Set Theme "Dracula"
Set TypingSpeed 35ms
Set Framerate 12

Hide
Type "export PATH=$PWD/build:$PATH && rm -rf /tmp/radar-demo && clear"
Enter
Show

Type "radar_sim info scenarios/default.ini"
Sleep 300ms
Enter
Sleep 3s

Type "clear && radar_sim run scenarios/default.ini -q -o /tmp/radar-demo"
Sleep 300ms
Enter
Sleep 5s

Type "clear && radar_sim montecarlo scenarios/monte_carlo.ini --trials 300 -o /tmp/radar-demo | head -21"
Sleep 300ms
Enter
Sleep 5s
