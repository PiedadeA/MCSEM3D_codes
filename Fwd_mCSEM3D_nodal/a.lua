-- Lua script.
p=tetview:new()
p:load_plc("modelo3d.poly")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
