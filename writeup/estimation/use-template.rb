#! /usr/bin/ruby

source = "estimation.tex"
destination = "estimation.tex"
template = "template.tex"

estimation = []
body = false

#Read the source file
File.open(source, "r").each_line do |line|
    if (!body and line.include?("\\section{Introduction}"))
        body = true
        puts "starting body"
    end

    if body
        if line.include?("\\bibliographystyle") 
            body = false
            puts "ending body"
        else
            estimation << line
        end
    end
end


#Read the template file:
template_contents = []
File.open(template, "r").each_line do |line|
    template_contents << line
end


#Replace the content section in template:
loc = template_contents.index {|s| s.include?("[content]")}
template_contents[loc..loc] = estimation


#Write the complete result to disk
File.open(destination, "w") do |f|
    template_contents.each {|element| f.write(element)}
end

