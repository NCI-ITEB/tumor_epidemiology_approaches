---
permalink: /instructors
layout: page
title: Instructors
nav_order: 4
---
<!----- <h1>Instructors</h1> ---->

  {% for author in site.authors %}
  <section>
	<figure style="float:left" class="image is-128x128">
		<img class="is-rounded" src="{{ site.baseurl }}/assets/img/{{ author.picture }}">
	</figure>
    <!---
    <img class="is-rounded" style="float:left" src="{{ site.baseurl }}/assets/img/{{ author.picture }}" alt="{{ author.name }}" height=180 width=180 />
    --->
    <div style="margin-left: 200px">
    	<h2>{{ author.name }}</h2>
    	<h4>{{ author.organization }}</h4>
    	<h4>{{ author.position }}</h4>
    	<p>{{ author.content | markdownify }}</p>
    </div>
  <br><br>
  {% endfor %}
  </section>
