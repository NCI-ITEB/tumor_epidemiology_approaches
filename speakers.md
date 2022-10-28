---
permalink: /speakers
layout: page
title: Invited Speakers
---

  <section class="px-4">
  {% for speaker in site.speakers %}
  <div>
    {% if speaker.picture %}
    <figure class="image is-128x128" style="float:left; clear:left" >
  		<img class="is-rounded" src="{{ site.baseurl }}/assets/speakers/{{ speaker.picture }}">
  	</figure>
    {% endif %}
    <div class="px-4" style="overflow:auto">
      	<h2>{{speaker.name}}, {{speaker.education}}</h2>
        <h5>{{speaker.affiliation}}</h5>
      	{% if speaker.talk_title %}
        <h4>{{ speaker.talk_date }} - <span class="is-italic">{{ speaker.talk_title }}</span></h4>
        {% else %}
        <h4>{{ speaker.talk_date }}</h4>
        {% endif %}
      	<p>{{ speaker.content | markdownify }}{% if speaker.profile_link %}<a href="{{speaker.profile_link}}" target="_blank">{{speaker.short_name}}'s profile</a> {% endif %}</p>
    </div>
    <br><br>
  </div>
  {% endfor %}
  </section>
