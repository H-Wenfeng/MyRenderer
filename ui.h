// #include <X11/Xlib.h>
// #include <X11/Xutil.h>

// void Create_Window()
// {
//     Display *display;
//     Window window;
//     XEvent event;
//     display = XOpenDisplay(NULL);
//     int s = DefaultScreen(display);
//     window = XCreateSimpleWindow(display, RootWindow(display, s), 10, 10, 900, 900, 1, BlackPixel(display, s), WhitePixel(display, s));
//     XSelectInput(display, window, ExposureMask | KeyPressMask);
//     XMapWindow(display, window);
//     GC gc = XCreateGC(display, window, 0, nullptr);
//     XFlush(display);
// };