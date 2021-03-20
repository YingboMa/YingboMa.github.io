<!--
Add here global page variables to use throughout your
website.
The website_* must be defined for the RSS to work
-->
@def website_title = "Yingbo Ma's homepage"
@def website_descr = "I put my random thoughts here."
@def website_url   = "https://YingboMa.github.io/"

@def author = "Yingbo Ma"

@def mintoclevel = 2

<!--
Add here files or directories that should be ignored by Franklin, otherwise
these files might be copied and, if markdown, processed by Franklin which
you might not want. Indicate directories by ending the name with a `/`.
-->
@def ignore = ["node_modules/", "franklin", "franklin.pub"]

<!--
Add here global latex commands to use throughout your
pages. It can be math commands but does not need to be.
For instance:
* \newcommand{\phrase}{This is a long phrase to copy.}
-->
\newcommand{\R}{\mathbb R}
\newcommand{\scal}[1]{\langle #1 \rangle}
\newcommand{\d}[1]{\mathop{}\!\mathrm{d} #1}
\newcommand{\dd}[2]{\frac{\d{#1}}{\d{#2}}}

\newcommand{\blogtitle}[1]{\#1}
\newcommand{\blogdate}[1]{~~~<p><span class="blog-date">#1</span></p>~~~}
