#! /usr/bin/ruby

#Parse command line arguments:
require 'optparse'

options = {}
OptionParser.new do |opts|
    opts.banner = "Usage: example.rb [options]"

    opts.on('-n', '--source SOURCE', 'Source file') { |v| options[:source] = v }
    opts.on('-h', '--destination DEST', 'Destination file') { |v| options[:destination] = v }
    opts.on('-p', '--template TEMPLATE', 'Template file') { |v| options[:template] = v }

end.parse!

estimation = []
proofs = []
body = false
appendix = false
theorem = false


#Read the source file
File.open(options[:source], "r").each_line do |line|
    if appendix
        if not line == "\n"
            if line.include?("\\bibliographystyle") 
                appendix = false
                puts "ending appendix"
            else
                proofs << line
            end
        end
    end
    
    #Dont include blank lines within theorems
    if body and theorem
        if line.include?("\\end{thm}")
            theorem = false
        elsif line != "\n"
            estimation << line
        end
    end

    if body and not theorem
        if line.include?("\\begin{thm}")
            theorem = true
            estimation << line
        elsif line.include?("\\appendix") 
            body = false
            appendix = true
            puts "ending body, beginning appendix"        
            proofs << line
        else
            estimation << line
        end
    end

    if (!body and !appendix and line.include?("\\begin{document}"))
        body = true
        puts "starting body"
    end
end


#Read the template file:
template_contents = []
File.open(options[:template], "r").each_line do |line|
    template_contents << line
end

#Replace the body section in template:
loc = template_contents.index {|s| s.include?("[body]")}
template_contents[loc..loc] = estimation

#Replace the appendix section in template:
loc = template_contents.index {|s| s.include?("[appendix]")}
template_contents[loc..loc] = proofs

#Write the complete result to disk
File.open(options[:destination], "w") do |f|
    template_contents.each {|element| f.write(element)}
end

