
#Uses the PIL/Pillow libraries
try:
    from PIL import Image, ImageDraw, ImageFont, ImageFilter
except:
    print("ERROR!  This script requires the Pillow library, which doesn't\
seem to be installed.")
    input("Press any key to quit.")
    exit()
    raise


#GFF file locations for both Windows and Linux (just to help development across machines)
linux_gff = [GFF LOCATION HERE]
windows_gff = [GFF LOCATION HERE]

#Query to search the GFF files for
query1 = "polyphos"
query2 = "ppk"


#Set colours for genes on positive and negative strands
positive_colour = (100,200,100)
negative_colour = (50,100,230)


def draw_arrow(image, xstart, ystart, length,strand,colour, linewidth=2, \
               arrow_width=25, arrowhead_length=25, ):

    #If the gene is on the negative strand flip the arrow
    if strand == "-":
        xstart = xstart + length
        length = -length
        arrowhead_length = -arrowhead_length
        

    #Shorten the gene length for the arrowhead (otherwise the arrow is drawn
    # as a rectangle the length of a gene plus an arrowhead)
    length = length - arrowhead_length
    
    #Draw the borderless triangle
    image.polygon([xstart + length, ystart - (arrow_width / 2),\
                  xstart + length + arrowhead_length, ystart + (arrow_width / 2),\
                  xstart + length, ystart + (arrow_width / 2) + arrow_width], \
                  fill = colour, outline=(0,0,0))
    
    #Draw the borderless rectangle
    image.rectangle([xstart, ystart, xstart + length, ystart + arrow_width],\
                   fill = colour)
    
    #Top horizontal outline
    image.line([xstart,ystart,xstart+length,ystart], fill=(0,0,0), \
              width=linewidth)
    #Top vertical outline
    image.line([xstart + length, ystart, xstart + length, ystart - (arrow_width / 2)], \
              fill=(0,0,0), width=linewidth)
    #Forward nose outline
    image.line([xstart + length, ystart - (arrow_width / 2), xstart + length + arrowhead_length,\
               ystart + (arrow_width / 2)], fill=(0,0,0), width=linewidth)
    #Reverse nose outline
    image.line([xstart + length + arrowhead_length, ystart + (arrow_width / 2), xstart + length,\
               ystart + (arrow_width / 2) + arrow_width],\
              fill=(0,0,0), width=linewidth)
    #Bottom vertical outline
    image.line([xstart + length, ystart + (arrow_width / 2) + arrow_width,\
               xstart + length, ystart + arrow_width], \
               fill=(0,0,0), width=linewidth)
    #Bottom horizontal outline
    image.line([xstart + length, ystart + arrow_width, xstart,\
               ystart + arrow_width], fill=(0,0,0),\
               width=linewidth)
    #Back vertical outline
    image.line([xstart, ystart + arrow_width, xstart, ystart],\
              fill=(0,0,0), width=linewidth)


print("Setting font...")
#Font locations for both Windows and Linux
#linux_font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/FreeSerif.ttf", 24)
windows_font = ImageFont.truetype("C:\\Windows\\winsxs\\arial_31bf3856ad364e35_6.1.7601.18528_none_d0a29012c3ff391b\\arial.ttf", 16)


print("Opening file...")
#Open the GFF and read it into the list called "lines"
with open(windows_gff, "r") as file:
    lines = file.readlines()

#List to hold data we're about to generate (dictionary of each seq in the GFF)
genes = []
#A counter to assign a number to each gene in the GFF.  This gives a numerical
#ID to each gene which makes it easier to index
counter = 0

print("Processing the GFF data...")
#For each gene in the input file:
for line in lines:
    #If the line is a gene, designated by "NODE" being on the line
    #if line[0:4] == "NODE":
    if (line[0] != "#" and line[0] != ">") and \
       (("Protein Homology" in line) or ("NODE" in line)):
        #print(repr(line))
    #if ("Protein Homology" in line) or ("NODE" in line) and ("#" not in line):
        #Split the line into a list of data (they're tab separated)
        split_data = line.split("\t")
        #print(split_data)
        #Create a dictionary describing the gene
        data_dictionary = {
            #First item is the counter
            "index":counter,
            "seqname":split_data[0],
            "source":split_data[1],
            "feature":split_data[2],
            "start":split_data[3],
            "end":split_data[4],
            "score":split_data[5],
            "strand":split_data[6],
            "frame":split_data[7],
            "attribute":split_data[8],
            }
        #Append the dictionary to the "data_lines" list
        genes.append(data_dictionary)
        #Increase the counter so that each sequence has a unique index
        counter += 1

indecies = []

print("Searching the GFF data for genes matching the query...")
#For each gene-dictionary
for seq in genes:
    #If the query is in the "attribute" data:
    if (query1.lower() in seq["attribute"].lower()) or \
       (query2.lower() in seq["attribute"].lower()):
        #record the index of the sequence so we can find it later
        indecies.append(seq["index"])

#print(indecies)
image_counter = 0

images = []
for index in indecies:
    images.append("img" + str(image_counter))
    image_counter = image_counter + 1

#print(images)

print("Drawing the diagrams...")
for position, index in enumerate(indecies):
    #A list to hold gene lengths
    lengths = []
    #A list to hold gene names
    names = []

    genes_of_interest = []

    #For 5 genes up and downstream of the gene of interest
    for x in range (index-5,  index+6):
        if genes[x]["seqname"] == genes[index]["seqname"]:
            genes_of_interest.append(genes[x])


    #Width of the image
    width = 1920
    #Height of the image
    height = 400

    #Get the total length of DNA to be represented, including gaps
    total_length = int(genes_of_interest[-1]["end"]) - int(genes_of_interest[0]["start"])

    #Need to scale the gene lengths for the width of the image, calculate this as
    #the percent of the total length of genes
    scale = total_length / width

    #For each gene, calculate the length and insert that into the gene dictionary
    for gene in genes_of_interest:
        length = (int(gene["end"]) - int(gene["start"])) / scale
        gene["length"] = length

    #For each gene, calculate the spacer and insert that into the gene dictionary
    for x in range (0, len(genes_of_interest)-1):
        spacer = (int(genes_of_interest[x+1]["start"]) - int(genes_of_interest[x]["end"])) / scale
        genes_of_interest[x]["spacer"] = spacer

    genes_of_interest[-1]["spacer"] = 0

    #Var to track how far over the image we've progressed with genes
    horiz_pos = 0

    #Var to track how far up text has gone
    vert_pos = 45
    drawing_plane = 350

    #Create an image with a white background
    images[position] = Image.new("RGB",(width,height), (255,255,255))
    #Draw the image.  Still don't understand why PIL works this way
    img2 = ImageDraw.Draw(images[position])

    #For each gene in "lengths"
    for x in range (0, len(genes_of_interest)):

        ###############
        #Add the label#
        ###############

        full_attrib = genes_of_interest[x]["attribute"]
        name_start = full_attrib.find("product=")
        name = full_attrib[name_start+8:]
        if ";" in name:
            name_end = name.find(";")
            name = name[:name_end]

        img2.text([horiz_pos, drawing_plane - vert_pos],name,\
                  font = windows_font,\
                  fill = (0,0,0))
        vert_pos = vert_pos + 30

        ################
        #Draw the arrow#
        ################

        if genes_of_interest[x]["strand"] == "+":
            colours = positive_colour
        elif genes_of_interest[x]["strand"] == "-":
            colours = negative_colour
        
        draw_arrow(img2, xstart=horiz_pos, ystart=drawing_plane,\
                   length=genes_of_interest[x]["length"],\
                   colour=colours,\
                   strand = genes_of_interest[x]["strand"])
        
        horiz_pos = (horiz_pos + genes_of_interest[x]["length"])

        ##################
        #Draw the spacers#
        ##################

        if x != (len(genes_of_interest) - 1):
            spacer = (int(genes_of_interest[x+1]["start"]) -\
                     int(genes_of_interest[x]["end"])) / scale
            img2.rectangle([horiz_pos, drawing_plane+11, horiz_pos + spacer,drawing_plane+13],\
                       fill = (0,0,0))
            horiz_pos = (horiz_pos + spacer)



        full_attrib = genes[index]["attribute"]
        name_start = full_attrib.find("product=")
        name = full_attrib[name_start+8:]
        if ";" in name:
            name_end = name.find(";")
            name = name[:name_end]

        title = windows_gff[::-1]
        end = title.find("\\")
        title  = title[:end]
        title = title[::-1]
        title = title + " " + name + " query= " + query1
        img2.text([5, 5],title,\
                  font = windows_font,\
                  fill = (0,0,0))


    #Show the image

print("Generating the composite image...")
number_images = len(images)
composite_height = (420 * number_images)

composite_image = Image.new("RGB",(1920,composite_height), (255,255,255))
for x in range (0,number_images):
    composite_image.paste(images[x],(0,420*x))

    
        
composite_image.show()
#img.save(composite_image + ".bmp", "bmp")
print("Done!")
