import riley
import kleinian
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import farey

###
# CONFIGURATION OPTIONS
###

limit_set_points = 2000
limit_set_depth = 5

scale = 100
riley_bounds = (-4,4,-4,4) # -x,x,-y,y
limit_bounds = (-4,4,-4,4) # -x,x,-y,y

slice_parabolic_max_denom = 30
slice_elliptic_max_denom = 30

## OPTIONS END


root = tk.Tk()
root.title("Riley slice limit sets")

def canvas_to_usual_coords(x,y):
    y = scale*(riley_bounds[3]-riley_bounds[2]) - y
    x = x/scale + riley_bounds[0]
    y = y/scale + riley_bounds[2]
    return x,y

def usual_coords_to_canvas(x,y):
    x = scale*(x - riley_bounds[0])
    y = scale*(y - riley_bounds[2])
    y = scale*(riley_bounds[3]-riley_bounds[2]) - y
    return x,y

selector_frame = ttk.Frame(root)
selector_frame.grid(column=0,row=0,sticky="nsew")
elliptic_p_var = tk.StringVar()
elliptic_p = 3
elliptic_p_var.set(str(elliptic_p))
elliptic_q_var = tk.StringVar()
elliptic_q = 4
elliptic_q_var.set(str(elliptic_q))
elliptic_frame = ttk.Frame(selector_frame)
elliptic_frame.grid(column=1,row=0)

p_entry=tk.Entry(elliptic_frame, textvariable=elliptic_p_var, state=tk.DISABLED)
q_entry=tk.Entry(elliptic_frame, textvariable=elliptic_q_var, state=tk.DISABLED)
slice_selection = tk.StringVar()
slice_selection.set('parabolic')
def change_slice(*args):
    new_slice = slice_selection.get()
    if new_slice == 'parabolic':
        redraw_slice(riley.riley_slice(np.inf,np.inf,slice_parabolic_max_denom))
        p_entry['state']=tk.DISABLED
        q_entry['state']=tk.DISABLED
    elif new_slice == 'elliptic':
        try:
            global elliptic_p, elliptic_q
            elliptic_p = int(elliptic_p_var.get())
            elliptic_q = int(elliptic_q_var.get())
            print('redraw with p = ' + str(elliptic_p) + ' q = ' + str(elliptic_q))
            redraw_slice(riley.riley_slice(elliptic_p,elliptic_q,slice_elliptic_max_denom))
            p_entry['state']=tk.NORMAL
            q_entry['state']=tk.NORMAL
        except ValueError:
            messagebox.showerror("Error", "Values for elliptic element orders must be integers")
    return True

p_entry.bind('<Return>', change_slice)
q_entry.bind('<Return>', change_slice)
tk.Radiobutton(selector_frame, text='Parabolic slice', variable=slice_selection, value='parabolic', command=change_slice).grid(column=0,row=0)
tk.Radiobutton(elliptic_frame, text='Elliptic slice', variable=slice_selection, value='elliptic', command=change_slice).pack()
tk.Label(elliptic_frame, text="p: ").pack(side='left')
p_entry.pack(side='left')
tk.Label(elliptic_frame, text="q: ").pack(side='left')
q_entry.pack(side='left')

mainframe = ttk.Frame(root)
mainframe.grid(column=0, row=1, sticky=(tk.N, tk.W, tk.E, tk.S))
slice_canvas = tk.Canvas(mainframe,width=scale*(riley_bounds[1]-riley_bounds[0]),height=scale*(riley_bounds[3]-riley_bounds[2]))
slice_canvas.grid(column=0,row=0)

def redraw_slice(points):
    slice_canvas.delete("all")
    for point in points:
        if point.real > riley_bounds[0]  and point.real < riley_bounds[1] and point.imag > riley_bounds[2] and point.imag < riley_bounds[3]:
            radius=1
            x,y = usual_coords_to_canvas(point.real,point.imag)
            slice_canvas.create_oval(x-radius, y-radius, x + radius, y + radius, fill="black", width=0)

change_slice()

positionframe = ttk.Frame(mainframe)
positionframe.grid(column=0, row=1, sticky='nsew')

hover_position = tk.StringVar()
ttk.Label(positionframe, textvariable=hover_position).pack()

selected_position = tk.StringVar()
ttk.Label(positionframe, textvariable=selected_position, foreground='red').pack()

limitset_canvas = tk.Canvas(mainframe,width=scale*(limit_bounds[1]-limit_bounds[0]),height=scale*(limit_bounds[3]-limit_bounds[2]))
limitset_canvas.grid(column=1,row=0)

mouse_down = False
def motion(event):
    x, y = event.x, event.y
    x,y = canvas_to_usual_coords(x,y)
    hover_position.set(str(x + y*1j))
    if mouse_down == True:
        print('redraw')
        redraw_limit(event.x,event.y)

last_selected = None
def redraw_limit(canvas_x,canvas_y):
    x,y = canvas_to_usual_coords(canvas_x,canvas_y)
    selected_position.set(str(x + y*1j))

    if auto_recompute.get() == 1:
        compute_farey()

    radius=5
    global last_selected
    if last_selected != None:
        slice_canvas.delete(last_selected)
    last_selected=slice_canvas.create_oval(canvas_x-radius,canvas_y-radius,canvas_x+radius,canvas_y+radius,fill='red',width=0)

    current_slice = slice_selection.get()
    if current_slice == 'parabolic':
        A = np.array([[1,1],[0,1]])
        B = np.array([[1,0],[x + y*1j,1]])
    elif current_slice == 'elliptic':
        A = np.array([[np.exp(2j*np.pi/elliptic_p),1],[0,np.exp(-2j*np.pi/elliptic_p)]])
        B = np.array([[np.exp(2j*np.pi/elliptic_q),0],[x + y*1j,np.exp(-2j*np.pi/elliptic_q)]])
    seed = farey.fixed_points(0,1,B[1][0],A[0][0],B[0][0])
    limitset = kleinian.limit_set_markov([A,B],seed,limit_set_depth,True,limit_set_points)
    colours = {-2: 'red', -1:'blue', 1:'green', 2:'yellow'}
    limitset_canvas.delete("all")
    for (points,colour) in limitset:
        for point in points:
          if point.real > limit_bounds[0]  and point.real < limit_bounds[1] and point.imag > limit_bounds[2] and point.imag < limit_bounds[3]:
              radius=1
              x,y = usual_coords_to_canvas(float(point.real),float(point.imag))
              #print(str(x) +' '+ str(y)  +' '+ str(point)+ ' '+ str(colour)+colours[colour])
              limitset_canvas.create_oval(x-radius, y-radius, x + radius, y + radius, fill=colours[colour], width = 0)

def mouse_click(event):
    global mouse_down
    mouse_down = True
    return redraw_limit(event.x,event.y)

def mouse_unclick(event):
    global mouse_down
    mouse_down = False

slice_canvas.bind('<ButtonPress>',mouse_click)
slice_canvas.bind('<ButtonRelease>',mouse_unclick)
slice_canvas.bind('<Motion>',motion)

def compute_farey(*args):
    try:
        vals = slope.get().split('/')
        if len(vals) != 2:
            raise ValueError
        word = farey.word(int(vals[0]),int(vals[1]))
        fareyword.set('W_'+slope.get()+' = ' + ''.join(word))

        if selected_position.get() != '':
            current_slice = slice_selection.get()
            if current_slice == 'parabolic':
                matrix = farey.matrix(int(vals[0]),int(vals[1]),complex(selected_position.get()),1,1)
            elif current_slice == 'elliptic':
                matrix = farey.matrix(int(vals[0]),int(vals[1]),complex(selected_position.get()),np.exp(2j*np.pi/elliptic_p),np.exp(2j*np.pi/elliptic_q))
            fareymatrix.set('matrix = ' + str(matrix))
            fareytrace.set('tr = ' + str(np.trace(matrix)))

    except ValueError:
        messagebox.showerror("Error", "Enter the slope in the format p/q with p,q integers")

farey_input_frame = ttk.Frame(mainframe)
farey_input_frame.grid(column=0,row=2)
ttk.Label(farey_input_frame, text="Slope in format p/q").grid(column=0, row=0, sticky=tk.W)
slope = tk.StringVar()
slope_entry = ttk.Entry(farey_input_frame, textvariable=slope, width=7)
slope_entry.grid(column=1, row=0, sticky=(tk.W,tk.E))
ttk.Button(farey_input_frame, text="Compute Farey word", command=compute_farey).grid(column=2, row=0, sticky=tk.W)

auto_recompute = tk.IntVar()
ttk.Checkbutton(farey_input_frame, text="Recompute on select?", variable=auto_recompute, onvalue=1, offvalue=0).grid(column=3, row=0, sticky=tk.E)

farey_output_frame = ttk.Frame(mainframe)
farey_output_frame.grid(column=1,row=2)
fareyword = tk.StringVar()
ttk.Label(farey_output_frame, textvariable=fareyword).grid(column=0, row=0, sticky=tk.W, padx=(0, 10))
fareymatrix = tk.StringVar()
ttk.Label(farey_output_frame, textvariable=fareymatrix).grid(column=1, row=0, padx=(0, 10), sticky=tk.W)
fareytrace = tk.StringVar()
ttk.Label(farey_output_frame, textvariable=fareytrace).grid(column=2, row=0, sticky=tk.E)



root.resizable(False, False)
root.mainloop()
