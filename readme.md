# Low-bandwidth video communication system
## Overview

This is the code from this [blog post](http://thomashansen.xyz/blog/low-bandwidth-video-streaming.html) 
which is looking into a method of sending video over a very small bandwidth 
(such as a few kilobytes). It accomplishes this by having pre-made sets of 
pixels on both the client and server, where each set has different concentrations
of pixels in different spots. The algorithm finds the best concentration of 
pixels to match the image, and then it only sends the pixels which are 
specified by the pre-made set. Once we've transmitted the pixels, we just 
interpolate the space inbetween the pixels, producing a lossy compression 
algorithm that is able to reduce images by a few 1000x.

Images tend to come out sofened, similar to running a median blur filter
over an image. The degree to which you reduce the image changes how soft the 
image appears.

Note: This is *not* a library, this is a proof of concept, which is an 
important distinction, because it means you will have to re write it to
implement it however you will be able to analyze how it runs easily from 
this code.

## Changes from the original write-up

This code has gone through a number of changes from the 
[original write-up](http://thomashansen.xyz/blog/low-bandwidth-video-streaming.html) 
I wrote about it. The main changes are how the best value_map is chose, where 
value_map referes to the pixels of each image that are sent (and specifically, 
which blocks are focused on when sending pixels). Now I use a rank,
where we rank the most important squares (high number = most important, low 
number = least important), and we compare those ranks to value maps using
immse. I'm still having a few issues with this method however I'm confident 
has the ability to be better than just using immse, especially for images that
have more rows and columns (i.e. higher resolution squares that partition the
image).

Note: I use the term "map" to refer to the set of pixels in an image, 
generally after a funciton has been applied to it (such as edge-detection, 
mean-square, or random generation).

Originally I had been trying to use immse or ssim functions to compare the 
different sets of pixels I could send, however I had issue with this since 
if the pixels in the pre-made map didn't line up correct frequently enough 
with the sobel pixels it wouldn't choose the map, which became an issue when 
you try to make the size of each seperate concentration smaller. The solution
to this is to rank each square (where each square has a different concentration),
and then compare the ranks, and not just the pixels.

```
ex. img           ex. closest-match set
 _ _ _ _            _ _ _ _ 
|5|8|3|4|          |4|7|2|5|
|6|1|2|7|          |6|1|3|8|

```
The ranking referes to the density of pixels in each square of the map, 
where a higher density earns a higher ranking. Currently each map is 
flattened and compared using immse(A,ref), with the closest match being used.

Additionally this is great because it allows you to use different methods to 
compare different maps. In this case you can now also use things like an LP 
since you're only comparing < 50 concentrations per image instead of millions
of pixels.
