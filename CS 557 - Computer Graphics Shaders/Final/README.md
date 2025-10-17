This GLSL fragment shader creates an animated celestial scene with fractal patterns and geometric shapes. Here's a technical breakdown of its core components and operations:

![alt text](https://github.com/Adoyle210/GradPortfolio/blob/main/CS%20557%20-%20Computer%20Graphics%20Shaders/Final/fractal3.png)

## Coordinate System Setup
uv = (gl_FragCoord.xy - 0.5 * u_resolution.xy) / min(u_resolution.x, u_resolution.y) * u_zoom;
I converted screen coordinates to normalized [-0.5,0.5] space and applied zoom control via u_zoom uniform. Use keys x and z to zoom in and out of the pattern. 
Signed Distance Functions (SDFs)
As mentioned in my project outline, I leverage iquilezles. For my fractal, I used the Moon, Star, and uneven capsule. 


## Fractal Pattern Generation
The shader employs a loop that iterates to generate intricate, self-similar patterns. This iterative process involves several key steps:
Coordinate fracturing is achieved by repeatedly scaling and wrapping the UV coordinates, creating a fragmented, repeating pattern across the screen.


An exponential falloff is applied based on the distance from the center, which adds depth and focus to the overall composition.


Sinusoidal modulation introduces a wave-like, time-dependent variation to the pattern, resulting in smooth, animated distortions.


The base fractal pattern is then artfully combined with the predefined celestial shapes (moon, star, and capsule) using a smooth minimum function, creating soft, organic transitions between different scene elements.


This multi-step process creates a complex, evolving visual tapestry that blends mathematical precision with artistic fluidity.
Color Dynamics
As mentioned in my project outline, I leveraged this article to create the color palette. The pallet generates color gradients using phase-shifted cosine waves. The parameters are as follows: Base offsets (a), Amplitude (b), frequency ( c), and phase shifts (d). 


## Rendering Pipeline
Distance field calculations for all shapes


Smooth blending between shape elements


Pattern sharpening via pow(0.01 / d, 1.2)


Color accumulation in loop iterations

## Questions
How does what you did compare with what you proposed?
	I followed and completed what I originally proposed as detailed above. 

What you didn't get working and why?
Using the guides/resources provided in my original outline I was able to get the full project to work. 

Anything extra you added?
I implemented key hits to zoom in and out of the fractal. (x/z keys) This allows for dynamic exploration of the fractal’s intricate details.  

Any particular cleverness you want me to know about?
The design lies in the hybrid rendering approach that merges the SDF primitives with fractal generation. Through strategic blending using smooth minimum operations, the system creates novel organic forms. 

What you learned by doing this?
	This assignment introduced me to the following: 
SDF construction and manipulation techniques
Fractal interaction through coordinate space transformations 
A crucial point was understanding how SDF combinations require careful distance field blending to maintain visual integrity. While initial shape merging produced unexpected results, this led to the discovery of more interesting visual phenomena through parameter experimentation. The process underscored how procedural systems can generate unanticipated complexity from simple mathematical operations. 

## Video
https://media.oregonstate.edu/media/t/1_09v3ksfz 
