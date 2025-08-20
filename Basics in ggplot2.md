# Basics-in-ggplot2

## Introduction

### Drawing your first plot
load the ggplot2 package
```
library(ggplot2)
```
Explore the mtcars data frame with str()
```
str(mtcars)
```
Execute the following command to create your first plot
```
ggplot(mrcars, aes(cyl, mpg))+
  geom_point()
```
Changing cyl to factor will change how the graph is display. Doing this will only present values that are present within the dataset.
```factor(cyl)```
### Mapping data columns to aesthetics
Edit to add a color aesthetic mapped to disp
```
ggplot(mtcars, aes(wt, mpg, color = disp)) +
  geom_point()
```
Try replacing color with size and see how that affects the plot that you just made. Assigning a variable to color or size will affect how the plot is displayed, assigning/grouping by color or size

### Adding geometrics
Explore the diamonds data frame with str()
```
str(diamonds)
```
Add ```geom_point()``` with ```+```. To add another line to a ggplot command a ```+``` has to go immediatly after the previous line is written, indicating that the command (or plot) is not done being created.
```
ggplot(diamonds, aes(carat, price)) +
  geom_point()
```
The ```geom_point()``` command will create a scatterplot. You can also add ```geom_smooth()``` to add a smooth trend line. Give it a shot by adding it to the plot that you just created.

#### Changing one geom or every geom
Map the color aesthetic to clarity
```
ggplot(diamonds, aes(carat, price, color = clarity)) +
  geom_point() +
  geom_smooth()
```
You can change the opacity of the points by adding ```alpha``` into the ```geom_point()``` command.
```
ggplot(diamonds, aes(carat, price, color = clarity)) +
  geom_point(alpha = 0.4) +
  geom_smooth()
```
This has change the opacity of the points by 40%.

#### Saving plots as vairables
When creating plots you are able to sae them as variables. By doing this, you can later on call the plot that you already made and add on additional characteristics later on. to do this, name you plot and follow it by a ```<-``` and the remaining code to create the plot.

## Aesthetics
The following link will provide additional information on how to change the shape of points and a link to how to gather specific color codes:

Point aesthetics: <https://r-charts.com/base-r/pch-symbols/#colors>
Color Codes: <https://www.rapidtables.com/web/color/RGB_Color.html>

### Visible Aesthetics

#### All about aesthetics: color, shape, and size
These are the aesthetics you can consider within aes() in this chapter: ```x```, ```y```, ```color```, ```fill```, ```size```, ```alpha```, ```labels``` and ```shape```.
#### All about aesthetics: color vs. fill
The ```color``` aesthetic changes the outline of a geom and the ```fill``` aesthetic changes the inside. ```geom_point()``` is an exception: you use ```color``` (not ```fill```) for the point color. However, some shapes have special behavior.
The default ```geom_point()``` uses ```shape = 19``` which is a solid circle. An alternative shape ```shape = 21```. This is a circle that allows you to use both ```fill``` for the inside and ```color``` for the outline. 
You can call ```?points()``` for a break down on the types of points.
```
ggplot(mtcars, aes(wt, mpg, fill = fcyl)) +
  geom_point(shape = 1, size = 4)
```
### Using attributes
#### All about attributes: color, shape, size, and alpha
You can specify colors in R using hex codes: a hash followed by two hexadecimal numbers each for red, green, and blue (```"#RRGGBB"```). Hexadecimal is base-16 counting. You have 0 to 9, and A representing 10 up to F representing 15. Pairs of hexadecimal numbers give you a range from 0 to 255. ```"#000000"``` is "black" (no color), ```"#FFFFFF"``` means "white", and ```"#00FFFF"``` is cyan (mixed green and blue).
### Modifying aesthetics
#### Updating aesthetic labels
There are two functions that can improve the look of your figure labels
```labs()``` will set the x- and y- axis labels and will take strings for each argument. ```scale_fill_manual()``` defines propers of the color scale (ie. axis). The first argument sets the legend title.

```
palette <- c(automatic = "#377EB8", manual = "#E41A1C")

#Set the position
ggplot(mtcars, aes(fcyl, fill = fam)) +
  geom_bar(position = "dodge") +
  labs(x = "Number of Cylinders", y = "Count")
  scale_fill_manual("Transmission", values = palette)
  ```
## Geometrics

### Scatter Plots

#### Overplotting 1: Large Datasets
Scatter plots can be created using 
```
# Plot price vs. carat, colored by clarity
plt_price_vs_carat_by_clarity <- ggplot(diamonds, aes(carat, price, color = clarity))

# Add a point layer with tiny points
plt_price_vs_carat_by_clarity + geom_point(alpha = 0.5, shape = ".")
```
#### Overplotting 2: Aligned Values
This is how to align values on a single axis. This occurs when one axis is continuous and the other is categorical, which can be overcome with some form of jittering 
```
# Plot base
plt_mpg_vs_fcyl_by_fam <- ggplot(mtcars, aes(fcyl, mpg, color = fam))

# Default points are shown for comparison
plt_mpg_vs_fcyl_by_fam + geom_point()
```
#### Overplotting 3: Low-precision data
You already saw how to deal with overplotting when using ```geom_point()``` in two cases:

Large datasets
Aligned values on a single axis
We used position = 'jitter' inside ```geom_point()``` or ```geom_jitter()```.

Let's take a look at another case:

Low-precision data
This results from low-resolution measurements like in the iris dataset, which is measured to 1mm precision (see viewer). It's similar to case 2, but in this case we can jitter on both the x and y axis.
```
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
  # Swap for jitter layer with width 0.1
  geom_jitter(alpha = 0.5, width = 0.1)
  
# To use a different approach
  # Set the position to jitter
  geom_point(alpha = 0.5, position = "jitter")
  
# To provide an alternative specification
 # Use a jitter position function with width 0.1
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1))
```
#### Overplotting 4: Integer data
This can be type ```integer``` (i.e. 1 ,2, 3…) or categorical (i.e. class ```factor```) variables. ```factor``` is just a special class of type ```integer```.

You'll typically have a small, defined number of intersections between two variables, which is similar to case 3, but you may miss it if you don't realize that integer and factor data are the same as low precision data.

The ```Vocab``` dataset provided contains the years of education and vocabulary test scores from respondents to US General Social Surveys from 1972-2004.
```
#Examine the structure of Vocab
str(Vocab)

#Plot vocabulary vs. education
ggplot(Vocab, aes(education, vocabulary)) +
  # Add a point layer
  geom_point()
  
 #Change to a jitter layer
  geom_jitter()
  
 #Set the transparency to 0.2
  geom_jitter(alpha = 0.2)

#Set the shape to 1
  geom_jitter(alpha = 0.2, shape = 1)
```
  
### Histograms

#### Drawing Histograms
Recall that histograms cut up a continuous variable into discrete bins and, by default, maps the internally calculated count variable (the number of observations in each bin) onto the y aesthetic. An internal variable called ```density``` can be accessed by using the ```..``` notation, i.e. ```..density...``` Plotting this variable will show the relative frequency, which is the height times the width of each bin.

```
#Plot mpg
ggplot(mtcars, aes(x=mpg)) +

  #Add a histogram layer
  geom_histogram()
  
   #Set the binwidth to 1
  geom_histogram(binwidth = 1)
  
  #Map y to ..density..
ggplot(mtcars, aes(mpg, ..density..)) +
  geom_histogram(binwidth = 1)
  datacamp_light_blue <- "#51A8C9"

ggplot(mtcars, aes(mpg, ..density..)) +
  #Set the fill color to datacamp_light_blue
  geom_histogram(binwidth = 1, fill = datacamp_light_blue)
```
#### Positions in histograms
Here, we'll examine the various ways of applying positions to histograms. ```geom_histogram()```, a special case of ```geom_bar()```, has a ```position``` argument that can take on the following values:

```stack``` (the default): Bars for different groups are stacked on top of each other.

```dodge```: Bars for different groups are placed side by side.

```fill```: Bars for different groups are shown as proportions.

```identity```: Plot the values as they appear in the dataset.

```
#Update the aesthetics so the fill color is by fam
ggplot(mtcars, aes(mpg, fill = fam)) +
  geom_histogram(binwidth = 1)
  
#Change the position to dodge
  geom_histogram(binwidth = 1, position = "dodge")
  
#Change the position to fill
  geom_histogram(binwidth = 1, position = "fill")
  
#Change the position to identity, with transparency 0.4
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.4)
```
### Bar Plots
#### Position in bar and col plots

Let's see how the position argument changes ```geom_bar()```.

We have three position options:

```stack```: The default
```dodge```: Preferred
```fill```: To show proportions

While we will be using ```geom_bar()``` here, note that the function ```geom_col()``` is just ```geom_bar()``` where both the position and stat arguments are set to ```"identity"```. It is used when we want the heights of the bars to represent the exact values in the data.

In this exercise, you'll draw the total count of cars having a given number of cylinders (```fcyl```), according to manual or automatic transmission type (```fam```).
```
#Plot fcyl, filled by fam
ggplot(mtcars, aes(fcyl, fill=fam)) +
  #Add a bar layer
  geom_bar()
  
  #Set the position to "fill"
  geom_bar(position = "fill")
  
  #Change the position to "dodge"
  geom_bar(position = "dodge")
````
#### Overlapping bar plots
You can customize bar plots further by adjusting the dodging so that your bars partially overlap each other. Instead of using position = "dodge", you're going to use position_dodge(), like you did with position_jitter() in the the previous exercises. Here, you'll save this as an object, posn_d, so that you can easily reuse it.

Remember, the reason you want to use position_dodge() (and position_jitter()) is to specify how much dodging (or jittering) you want.

For this example, you'll use the mtcars dataset.
```
ggplot(mtcars, aes(cyl, fill = fam)) +
  # Change position to use the functional form, with width 0.2
  geom_bar(position = position_dodge(width=0.2))

# Set the transparency to 0.6
  geom_bar(position = position_dodge(width = 0.2), alpha = 0.6)
```
#### Bar plots: Sequential color palette
In this bar plot, we'll fill each segment according to an ordinal variable. The best way to do that is with a sequential color palette.

Here's an example of using a sequential color palette with the mtcars dataset:

ggplot(mtcars, aes(fcyl, fill = fam)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set1")
In the exercise, you'll use similar code on the the Vocab dataset. Both datasets are ordinal.

```
# Plot education, filled by vocabulary
ggplot(Vocab, aes(education, fill = vocabulary)) +
  # Add a bar layer with position "fill"
  geom_bar(position = "fill") +
  # Add a brewer fill scale with default palette
  scale_fill_brewer()
```
### Line Plots
#### Basic Line Plots
```
# Print the head of economics
head(economics)

# Using economics, plot unemploy vs. date
ggplot(economics, aes(date, unemploy)) +
  # Make it a line plot

# Change the y-axis to the proportion of the population that is unemployed
ggplot(economics, aes(date, unemploy/pop)) +
  geom_line()
```
#### Multiple time Series
```
# Plot the Rainbow Salmon time series
ggplot(fish.species, aes(x = Year, y = Rainbow)) +
  geom_line()

# Plot the Pink Salmon time series
ggplot(fish.species, aes(x = Year, y = Pink)) +
  geom_line()

# Plot multiple time-series by grouping by species
ggplot(fish.tidy, aes(Year, Capture)) +
  geom_line(aes(group = Species))+
```
## Themes
### Themes From Scratch
#### Moving the legend
Let's wrap up this course by making a publication-ready plot communicating a clear message.

To change stylistic elements of a plot, call theme() and set plot properties to a new value. For example, the following changes the legend position.

p + theme(legend.position = new_value)
Here, the new value can be

"top", "bottom", "left", or "right'": place it at that side of the plot.
"none": don't draw it.
c(x, y): c(0, 0) means the bottom-left and c(1, 1) means the top-right.
Let's revisit the recession period line plot (assigned to plt_prop_unemployed_over_time).

```
# View the default plot
plt_prop_unemployed_over_time

# Remove legend entirely
plt_prop_unemployed_over_time +
  theme(legend.position = "none")
  
# Position the legend at the bottom of the plot
plt_prop_unemployed_over_time +
  theme(legend.position = "bottom")
  
# Position the legend inside the plot at (0.6, 0.1)
plt_prop_unemployed_over_time +
  theme(legend.position = c(0.6,0.1))
```

#### Modifying theme elements
Many plot elements have multiple properties that can be set. For example, line elements in the plot such as axes and gridlines have a color, a thickness (size), and a line type (solid line, dashed, or dotted). To set the style of a line, you use element_line(). For example, to make the axis lines into red, dashed lines, you would use the following.

p + theme(axis.line = element_line(color = "red", linetype = "dashed"))
Similarly, element_rect() changes rectangles and element_text() changes text. You can remove a plot element using element_blank().

plt_prop_unemployed_over_time is available.

```
plt_prop_unemployed_over_time +
  theme(
    rect = element_rect(fill = "grey92"),
    legend.key = element_rect(color = NA),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    # Add major y-axis panel grid lines back
    panel.grid.major.y = element_line(
      # Set the color to white
      color = "white",
      # Set the size to 0.5
      size = 0.5,
      # Set the line type to dotted
      linetype = "dotted"
    )
  )
```

#### Modifying whitespace
Whitespace means all the non-visible margins and spacing in the plot.

To set a single whitespace value, use unit(x, unit), where x is the amount and unit is the unit of measure.

Borders require you to set 4 positions, so use margin(top, right, bottom, left, unit). To remember the margin order, think TRouBLe.

The default unit is "pt" (points), which scales well with text. Other options include "cm", "in" (inches) and "lines" (of text).

plt_mpg_vs_wt_by_cyl is available. The panel and legend are wrapped in blue boxes so you can see how they change.

```
# View the original plot
plt_mpg_vs_wt_by_cyl

plt_mpg_vs_wt_by_cyl +
  theme(
    # Set the axis tick length to 2 lines
    axis.ticks.length = unit(2, "lines")
  )
   # Set the legend key size to 3 centimeters
    legend.key.size = unit(3, "cm")
  )
# Set the legend margin to (20, 30, 40, 50) points
    legend.margin = unit(margin(20, 30, 40, 50, "pt")
  ))
# Set the plot margin to (10, 30, 50, 70) millimeters
    plot.margin =unit(margin(10,30,50,70,"mm")
  ))
```

### Theme Flexibility
#### Built-in Themes
In addition to making your own themes, there are several out-of-the-box solutions that may save you lots of time.

theme_gray() is the default.
theme_bw() is useful when you use transparency.
theme_classic() is more traditional.
theme_void() removes everything but the data.
plt_prop_unemployed_over_time is available.

```
# Add a black and white theme
plt_prop_unemployed_over_time +
  theme_bw()

# Add a classic theme
plt_prop_unemployed_over_time +
  theme_classic()
  
# Add a void theme
plt_prop_unemployed_over_time +
  theme_void()
```

#### Exploring ggthemes
Outside of ggplot2, another source of built-in themes is the ggthemes package. The workspace already contains the plt_prop_unemployed_over_time, the line plot from before. Let's explore some of the ready-made ggthemes themes.

plt_prop_unemployed_over_time is available.

```
# Use the fivethirtyeight theme
plt_prop_unemployed_over_time +
  theme_fivethirtyeight()

# Use Tufte's theme
plt_prop_unemployed_over_time +
  theme_tufte()

# Use the Wall Street Journal theme
plt_prop_unemployed_over_time +
  theme_wsj()
```

#### Setting themes
Reusing a theme across many plots helps to provide a consistent style. You have several options for this.

Assign the theme to a variable, and add it to each plot.
Set your theme as the default using theme_set().
A good strategy that you'll use here is to begin with a built-in theme then modify it.

plt_prop_unemployed_over_time is available. The theme you made earlier is shown in the sample code.

```
theme_recession <- theme(
  rect = element_rect(fill = "grey92"),
  legend.key = element_rect(color = NA),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.grid.major.y = element_line(color = "white", size = 0.5, linetype = "dotted"),
  axis.text = element_text(color = "grey25"),
  plot.title = element_text(face = "italic", size = 16),
  legend.position = c(0.6, 0.1)
)
theme_tufte_recession <- theme_tufte() + theme_recession

#Set theme_tufte_recession as the default theme
theme_set(theme_tufte_recession)

#Draw the plot (without explicitly adding a theme)
plt_prop_unemployed_over_time
```

#### Publication-quality plots
```
plt_prop_unemployed_over_time +
  theme_tufte() +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    # Set the axis title's text color to grey60
    axis.title = element_text(color = "grey60"),
    # Set the axis text's text color to grey60
    axis.text = element_text(color = "grey60")
  )
```
### Effective explanatory plots
#### Using geoms for explanatory plots
Let's focus on producing beautiful and effective explanatory plots. In the next couple of exercises, you'll create a plot that is similar to the one shown in the video using gm2007, a filtered subset of the gapminder dataset.

This type of plot will be in an info-viz style, meaning that it would be similar to something you'd see in a magazine or website for a mostly lay audience.

A scatterplot of lifeExp by country, colored by lifeExp, with points of size 4, is provided.

```
#Set the color scale
palette <- brewer.pal(5, "RdYlBu")[-(2:4)]

#Modify the scales
ggplot(gm2007, aes(x = lifeExp, y = country, color = lifeExp)) +
  geom_point(size = 4) +
  geom_segment(aes(xend = 30, yend = country), size = 2) +
  geom_text(aes(label = round(lifeExp,1)), color = "white", size = 1.5) +
  scale_x_continuous("", expand = c(0,0), limits = c(30,90), position = "top") +
  scale_color_gradientn(colors = palette)
```
#### Using annotate() for embellishments
In the previous exercise, we completed our basic plot. Now let's polish it by playing with the theme and adding annotations. In this exercise, you'll use ```annotate()``` to add text and a curve to the plot.

The following values have been calculated for you to assist with adding embellishments to the plot:

```
global_mean <- mean(gm2007_full$lifeExp)
x_start <- global_mean + 4
y_start <- 5.5
x_end <- global_mean
y_end <- 7.5
Our previous plot has been assigned to plt_country_vs_lifeExp
```

```
#Define the theme
plt_country_vs_lifeExp +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        legend.position = "none")
        
#Add a vertical line
plt_country_vs_lifeExp +
  step_1_themes +
  geom_vline(xintercept = global_mean, color = "grey40", linetype = 3)

#Add text
plt_country_vs_lifeExp +
  step_1_themes +
  geom_vline(xintercept = global_mean, color = "grey40", linetype = 3) +
  annotate(
    "text",
    x = x_start, y = y_start,
    label = "The\nglobal\naverage",
    vjust = 1, size = 3, color = "grey40"
  )
```
