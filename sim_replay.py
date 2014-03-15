# imports -----------------------------------------
import OpenGL
from OpenGL.GLU import gluLookAt
OpenGL.ERROR_ON_COPY = True
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.constants import GLfloat

import sys, os, random, time
from math import sin, cos, sqrt, pi

import numpy as np

# view settings
(view_rotx, view_roty, view_dist)=(-40.0, 90.0, 2000.0)

# which index from data to draw
curIndex = 0
# which datafile to draw from
curData = 0


# change view angle, exit upon ESC
def key(k, x, y):
    global view_dist, curIndex, curData
    
    if k == 'a':
        view_dist *= 0.9
    elif k == 'z':
        view_dist *= 1.1
    elif k == ' ':
        curIndex += 1
    elif k in [str(i) for i in range(1,nfiles+1)]:
        curData = int(k)-1
    elif ord(k) == 27: # Escape
        glutDestroyWindow(glutGetWindow())
        glutLeaveMainLoop()
        return
    else:
        return
    glutPostRedisplay()


# change view angle
def special(k, x, y):
    global view_rotx, view_roty

    if k == GLUT_KEY_UP:
        view_rotx += 5.0
    elif k == GLUT_KEY_DOWN:
        view_rotx -= 5.0
    elif k == GLUT_KEY_LEFT:
        view_roty += 5.0
    elif k == GLUT_KEY_RIGHT:
        view_roty -= 5.0
    else:
        return
    glutPostRedisplay()


# new window size or exposure
def reshape(width, height):
    h = float(height) / float(width);
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glFrustum(-1.0, 1.0, -h, h, 50.0, 50000.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

def initGL():
    # set various openGL params
    glDisable(GL_CULL_FACE)
    glDisable(GL_DEPTH_TEST)

    glEnable (GL_BLEND)
    glBlendFunc (GL_SRC_ALPHA, GL_ONE)

    glEnable(GL_NORMALIZE)


def draw():
    global curIndex, curData, data
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    if curData>2 or curData==0:
        pz = 400
    else:
        pz = 0
    glPushMatrix()
    gluLookAt(0, view_dist, 0, 0, 0, 0, 0, 0, 1);
    glRotatef(view_rotx, 1, 0, 0)
    glRotatef(view_roty, 0, 0, 1)
    glTranslatef(0,0,-pz)
    
    # then the DNA
    glColor(1,1,1,1)
    glBegin(GL_LINE_STRIP)
    if curIndex >= data[curData].shape[0]:
        curIndex = 0
    for vert in data[curData][curIndex,:,:]:
        glVertex3f(vert[0],vert[1],vert[2])
    glEnd()
    
    # and a pore, if needed
    if curData > 4:
        glColor(1,0,0,1)
        r = 5
        z = 400
        
        glBegin(GL_LINE_STRIP)
        for th in np.linspace(0,2*np.pi,20):
            glVertex3f(r*np.cos(th),r*np.sin(th),z)
        glEnd()
    
    glPopMatrix()

    # now draw things in screen coords
    # GL projects rasterpos, so need to set appropriate proj mat
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0,1,0,1,-1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glColor(1,1,1,1)
    glRasterPos2f(0.02, 0.02);
    # loop through string and spit out per char
    debugstr = "Index: " + str(curIndex)
    for c in debugstr:
       glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(c));

    glMatrixMode(GL_PROJECTION);
    glPopMatrix()
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix()

    glutSwapBuffers()
    
    glutPostRedisplay()


def gogogo():
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)

    glutInitWindowPosition(0, 0)
    glutInitWindowSize(500, 500)
    glutCreateWindow("DNA Monte Carlo")
    
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION) 
    glutSetOption(GLUT_ACTION_GLUTMAINLOOP_RETURNS, GLUT_ACTION_CONTINUE_EXECUTION) 

    initGL()

    glutDisplayFunc(draw)
    glutReshapeFunc(reshape)
    glutKeyboardFunc(key)
    glutSpecialFunc(special)
    
    glutMainLoop()


# load files
nfiles = 7
data = []
for i in range(nfiles):
    dat=np.load('data/run'+str(i)+'.npy')
    data.append(dat)
    
gogogo()